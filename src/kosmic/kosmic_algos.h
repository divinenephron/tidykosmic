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

namespace kosmic {

	/*
	 Calculates and evaluates estimated normal distributions.
	 Requires an internal est_cdf array for temporary calculation.
	*/
	template<typename T> class abstract_cost_evaluator {
	protected:

		/*
		 Evaluates the given truncation interval t1 <= x < t2 and the intervals < and >= the truncation interval.
		*/
		virtual T cost(int t1_i, int t2_i) const = 0;

		/*
		 Evaluates all possible (with respect to t1/t2 min/max) truncation intervals and returns the interval with the
		 minimum penalty.
		*/
		evaluated_interval<T> min_cost() const {
			evaluated_interval<T> best_cost(INFINITY);
			for(int t1_i = t1min_i; t1_i <= t1max_i; t1_i++) {
				for(int t2_i = t2min_i; t2_i <= t2max_i; t2_i++) {
					T temp_cost = cost(t1_i, t2_i);
					if(temp_cost < best_cost.cost)
						best_cost = evaluated_interval<T>(temp_cost, t1_i, t2_i);
				}
			}
			return best_cost;
		}

		int classes;
		const T* cdf_emp;
		T* cdf_est;
		const T* x;
		int t1min_i, t1max_i, t2min_i, t2max_i;

	public:

		abstract_cost_evaluator(int classes, const T cdf_emp[], T cdf_est[], const T x[], int t1min_i, int t1max_i, int t2min_i, int t2max_i) :
			classes(classes), cdf_emp(cdf_emp), cdf_est(cdf_est), x(x), t1min_i(t1min_i), t1max_i(t1max_i), t2min_i(t2min_i), t2max_i(t2max_i) {
		}

		/*
		 Calculates and evaluates cdf_est for mu/sigma.
		*/
		evaluated_interval<T> evaluate(T mu, T sigma) {
			if(sigma <= 0)
				return INFINITY;
			for(int x_i = 0; x_i < classes; x_i++)
				cdf_est[x_i] = std_normal_cdf<T>(mu, sigma, x[x_i]);
			return min_cost();
		}

	};

	/*
	 Standard cost function, as per https://www.nature.com/articles/s41598-020-58749-2.
	*/
	template<typename T> class cost_evaluator_ks_classic : public abstract_cost_evaluator<T> {

		T cost(int t1_i, int t2_i) const override {
			auto F_emp = [this](int x_i) -> T { return this->cdf_emp[x_i]; };
			T a = 1 / (this->cdf_est[t2_i] - this->cdf_est[t1_i]) * (F_emp(t2_i) - F_emp(t1_i));
			T b = -this->cdf_est[t1_i] * a + F_emp(t1_i);
			auto F_est = [this, a, b](int x_i) -> T { return this->cdf_est[x_i] * a + b; };
			T max_ks_1 = 0.0;
			for(int i = 0; i < t1_i; i++)
				max_ks_1 = std::max(max_ks_1, F_emp(i) - F_est(i));
			T max_ks_t = 0.0;
			for(int i = t1_i; i < t2_i; i++)
				max_ks_t = std::max(max_ks_t, std::abs(F_est(i) - F_emp(i)));
			T max_ks_2 = 0.0;
			for(int i = t2_i; i < this->classes; i++)
				max_ks_2 = std::max(max_ks_2, F_est(i) - F_emp(i));
			return (max_ks_t + max_ks_1 + max_ks_2) / std::sqrt(F_est(t2_i) - F_est(t1_i));
		}

	public:
		cost_evaluator_ks_classic(int classes, const T cdf_emp[], T cdf_est[], const T x[], int t1min_i, int t1max_i, int t2min_i, int t2max_i) :
			abstract_cost_evaluator<T>(classes, cdf_emp, cdf_est, x, t1min_i, t1max_i, t2min_i, t2max_i) {}

	};

	/*
	 Improved cost function.
	*/
	template<typename T> class cost_evaluator_ks_log : public abstract_cost_evaluator<T> {

		T cost(int t1_i, int t2_i) const override {
			auto F_emp = [this](int x_i) -> T { return this->cdf_emp[x_i]; };
			T a = 1 / (this->cdf_est[t2_i] - this->cdf_est[t1_i]) * (F_emp(t2_i) - F_emp(t1_i));
			T b = -this->cdf_est[t1_i] * a + F_emp(t1_i);
			auto F_est = [this, a, b](int x_i) -> T { return this->cdf_est[x_i] * a + b; };
			T max_ks_1 = 0.0;
			for(int i = 0; i < t1_i; i++)
				max_ks_1 = std::max(max_ks_1, F_emp(i) - F_est(i));
			T max_ks_t = 0.0;
			for(int i = t1_i; i < t2_i; i++)
				max_ks_t = std::max(max_ks_t, std::abs(F_est(i) - F_emp(i)));
			T max_ks_2 = 0.0;
			for(int i = t2_i; i < this->classes; i++)
				max_ks_2 = std::max(max_ks_2, F_est(i) - F_emp(i));
			return (max_ks_t + max_ks_1 + max_ks_2) / std::log(1 + F_est(t2_i) - F_est(t1_i));
		}

	public:
		cost_evaluator_ks_log(int classes, const T cdf_emp[], T cdf_est[], const T x[], int t1min_i, int t1max_i, int t2min_i, int t2max_i) :
			abstract_cost_evaluator<T>(classes, cdf_emp, cdf_est, x, t1min_i, t1max_i, t2min_i, t2max_i) {}

	};

}


