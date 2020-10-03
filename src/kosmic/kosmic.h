/*
 Reference Interval Estimation from Mixed Distributions using Truncation Points
 and the Kolmogorov-Smirnov Distance (kosmic)

 Copyright (c) 2020 Jakob Zierk
 jakobtobiaszierk@gmail.com or jakob.zierk@uk-erlangen.de

 Based on ideas and prior implementations by Farhad Arzideh
 farhad.arzideh@gmail.com or farhad.arzideh@uk-koeln.de

 This program is free software: you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation, either version 3 of the License, or (at your option) any later
 version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with
 this program. If not, see <https://www.gnu.org/licenses/>.
*/
#pragma once

#include <algorithm>
#include <cmath>

#include "../lib/ctpl_stl.h"
#include "../hist_sampler_r.h"

#include "kosmic_structs.h"
#include "kosmic_math.h"
#include "kosmic_algos.h"
#include "kosmic_input.h"
#include "nm2_optimizer.h"

namespace kosmic {

	/*
	 kosmic main class (reference interval estimation algorithm).
	 Estimation of reference limits using truncation.
	 The truncation interval is defined as x >= t1, x < t2.
	 t1/t2 are optimized within t1/t2 >= t1min/t2min, t1/t2 <= t1max/t2max.
	*/
	template<typename T, typename cost_evaluator> class ri_estimator {

		struct {
			int classes;
			T* cdf_emp;
			T* x;
			int t1min_i, t1max_i, t2min_i, t2max_i;
			int mode_i, sd_i;
			double abs_tol;
		} data;

	public:

		/*
		 Initializes kosmic.
		*/
		ri_estimator(int classes, T* x, T* cdf, const algorithm_settings& settings) {
			data.classes = classes;
			data.x = x;
			data.cdf_emp = cdf;

			// Find percentiles according to settings:
			// Possible optimization: replace w/ lower_bound etc.
			data.t1min_i = data.t1max_i = data.t2min_i = data.t2max_i = -1;
			data.sd_i = -1;
			for(int i = 0; i < data.classes - 1; i++) {
				if(data.t1min_i == -1 && data.cdf_emp[i + 1] >= settings.t1min_p)
					data.t1min_i = i;
				if(data.t1max_i == -1 && data.cdf_emp[i + 1] > settings.t1max_p)
					data.t1max_i = i;
				if(data.sd_i == -1 && data.cdf_emp[i] >= settings.sd_p)
					data.sd_i = i;
				if(data.t2min_i == -1 && data.cdf_emp[i + 1] >= settings.t2min_p)
					data.t2min_i = i;
				if(data.t2max_i == -1 && data.cdf_emp[i + 1] > settings.t2max_p)
					data.t2max_i = i;
			}

			// Find mode:
			T max_ratio = data.cdf_emp[0];
			for(int i = 1; i < data.classes - 1; i++) {
				if(data.cdf_emp[i + 1] - data.cdf_emp[i] > max_ratio) {
					max_ratio = data.cdf_emp[i + 1] - data.cdf_emp[i];
					data.mode_i = i + 1;
				}
			}

			data.abs_tol = settings.abs_tol;
		}

		/*
		 Checks whether the input dataset is valid.
		*/
		bool is_good() {
			if(!(
				data.classes >= 6 && data.t2max_i != -1 && data.sd_i != -1 &&
				data.t1min_i < data.t1max_i && data.t1max_i < data.t2min_i - 1 && data.t2min_i < data.t2max_i &&
				data.mode_i != 0 && data.mode_i != data.classes - 1
				))
				return false;
			return true;
		}

		/*
		 Determines the best Normal distribution for data.x.
		 Overwrites the input arrays x[] and cdf_est of classes data.classes.
		*/
		best_normal<T> find_best_normal(T cdf_est[], T x[]) {
			cost_evaluator e(data.classes, data.cdf_emp, cdf_est, x, data.t1min_i, data.t1max_i, data.t2min_i, data.t2max_i);
			auto opt = nm2_optimizer<T, evaluated_interval<T>, cost_evaluator>::optimize(
				x[data.mode_i], x[data.sd_i] - x[data.mode_i],
				.95 * x[data.mode_i], x[data.mode_i + 1] - x[data.mode_i - 1],
				1.05 * x[data.mode_i], (x[data.mode_i + 1] - x[data.mode_i - 1] + x[data.sd_i] - x[data.mode_i]) / 2,
				e, data.abs_tol);
			return best_normal<T>{ opt.best_x0, opt.best_x1, opt.best_objective };
		}

		/*
		 Convenience method for find_best_normal(T x[], T cdf_est[]), creates and frees temporary arrays x[] and cdf_est[].
		*/
		best_normal<T> find_best_normal() {
			T* cdf_est = new T[data.n];
			best_normal<T> res = find_best_normal(cdf_est, data.x);
			delete[] cdf_est;
			return res;
		}

		/*
		 Determines the best Normal distribution after Box Cox transformation of data.x.
		 Uses Brute Force to determine the optimal lambda (0.0, 0.1, ... 1.0, and
		 optimum-0.09, optimum-0.08, ... optimum+0.09 afterwards).
		 Overwrites the input arrays x[] and cdf_est[] (length = data.classes).
		*/
		best_box_cox_normal<T> find_best_normal_after_boxcox(T x[], T cdf_est[]) {
			best_normal<T> best;
			best.interval.cost = INFINITY;
			T best_lambda;
			for(T lambda = 0.0; lambda <= 1.0; lambda += 0.1) {
				for(int i = 0; i < data.classes; i++)
					x[i] = box_cox_transform<T>(data.x[i], lambda);
				best_normal<T> bn = find_best_normal(cdf_est, x);
				if(bn.interval.cost < best.interval.cost) {
					best = bn;
					best_lambda = lambda;
				}
			}
			T best_lambda_2 = best_lambda;
			for(T lambda = std::max(0.0, best_lambda - 0.09); lambda < std::min(1.0, best_lambda + 0.1); lambda += 0.01) {
				if(lambda != best_lambda) {
					for(int i = 0; i < data.classes; i++)
						x[i] = box_cox_transform<T>(data.x[i], lambda);
					best_normal<T> bn = find_best_normal(cdf_est, x);
					if(bn.interval.cost < best.interval.cost) {
						best = bn;
						best_lambda_2 = lambda;
					}
				}
			}
			return best_box_cox_normal<T>{ best_lambda_2, best.mu, best.sigma, best.interval };
		}

		/*
		 Convenience method for find_best_normal_after_boxcox(T x[], T cdf_est[]), creates and frees temporary arrays x[] and cdf_est[].
		*/
		best_box_cox_normal<T> find_best_normal_after_boxcox() {
			T* x = new T[data.classes];
			T* cdf_est = new T[data.classes];
			best_box_cox_normal<T> res = find_best_normal_after_boxcox(x, cdf_est);
			delete[] x;
			delete[] cdf_est;
			return res;
		}

		/*
		 Run kosmic algorithm.
		 Return:
		  0: Success.
		  1: kosmic::is_good() == false
		  2: any exception
		*/
		static int run(const hist_builder<T> & hist, const algorithm_settings & settings,
			best_box_cox_normal<T> & result, best_box_cox_normal<T> * bootstrap_results, int bootstrap, int bootstrap_seed,
			int threads) {
			assert(bootstrap >= 0);
			assert(bootstrap_seed >= 0);
			assert(threads > 0);

			try {
				cdf<T> cdf(hist);

				// Multithread or single-thread version:
				if(threads > 1 && bootstrap > 0) {
					ctpl::thread_pool pool(threads);

					// Main result RI estimation:
					bool data_is_good = true;
					pool.push([&](int id) {
						ri_estimator<T, cost_evaluator> ri_estimator(cdf.classes, cdf.x, cdf.cum_freq, settings);
						if(!ri_estimator.is_good())
							data_is_good = false;
						else
							result = ri_estimator.find_best_normal_after_boxcox();
						});

					// RI bootstrapping:
					hist_sampler_r<T> hist_sampler(hist);
					int* temp_counts = new int[threads * hist.classes()];
					T* temp_cdf_emps = new T[threads * hist.classes()];
					T* temp_xs = new T[threads * hist.classes()];
					T* temp_cdf_ests = new T[threads * hist.classes()];
					for(int i = 0; i < bootstrap; i++) {
						pool.push([=, &hist_sampler, &cdf, &settings, &data_is_good](int id) {
							// Create random samples, might fail (bad sample due to random effects -> try again, or bad dataset -> abort):
							for(int j = 0;; j += bootstrap) {
								hist_sampler.cdf(i + j + bootstrap_seed, &temp_counts[id * hist_sampler.classes()], &temp_cdf_emps[id * hist_sampler.classes()]);
								ri_estimator<T, cost_evaluator> ri_estimator(hist_sampler.classes(), cdf.x, &temp_cdf_emps[id * hist_sampler.classes()], settings);
								if(!ri_estimator.is_good()) {
									if(!data_is_good)
										return;
									else
										continue;
								}
								bootstrap_results[i] = ri_estimator.find_best_normal_after_boxcox(&temp_xs[id * hist_sampler.classes()], &temp_cdf_ests[id * hist_sampler.classes()]);
								break;
							}
							});
					}

					// Wait for all threads to execute and delete temporary memory:
					pool.stop(true);
					delete[] temp_counts;
					delete[] temp_cdf_emps;
					delete[] temp_xs;
					delete[] temp_cdf_ests;

					return data_is_good ? 0 : 1;
				} else {
					ri_estimator<T, cost_evaluator> ri_estimator_(cdf.classes, cdf.x, cdf.cum_freq, settings);
					if(!ri_estimator_.is_good())
						return 1;

					// Main result RI estimation:
					result = ri_estimator_.find_best_normal_after_boxcox();

					// RI bootstrapping:
					if(bootstrap > 0) {
						int* temp_counts = new int[hist.classes()];
						T* temp_cdfs = new T[hist.classes()];
						hist_sampler_r<T> hist_sampler(hist);
						for(int i = 0; i < bootstrap; i++) {
							hist_sampler.cdf(i, temp_counts, temp_cdfs);
							// Create random samples, might fail (bad sample due to random effects -> try again):
							for(;;) {
								ri_estimator<T, cost_evaluator> ri_estimator_(hist_sampler.classes(), cdf.x, temp_cdfs, settings);
								if(!ri_estimator_.is_good())
									continue;
								bootstrap_results[i] = ri_estimator_.find_best_normal_after_boxcox();
								break;
							}
						}
						delete[] temp_counts;
						delete[] temp_cdfs;
					}

					return 0;
				}
			} catch(...) {
				return 2;
			}
		}

	};

}

