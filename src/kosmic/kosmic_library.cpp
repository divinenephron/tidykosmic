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

#ifdef _WINDOWS
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif

#include "kosmic.h"
#include "kosmic_algos.h"
#include "kosmic_library.h"
#include "kosmic_input.h"

using namespace kosmic;

#ifdef _WINDOWS
BOOL APIENTRY DllMain(HMODULE hModule, DWORD ul_reason_for_call, LPVOID lpReserved) {
	return TRUE;
}
#endif

#ifdef _WINDOWS
extern "C" __declspec(dllexport)
#else
extern "C"
#endif
int kosmic_lib(const double input[], int n, int decimals, int bootstrap, int bootstrap_seed, int threads, 
	double t1min, double t1max, double t2min, double t2max, double sd, double tol, char* _cost_func,
	best_box_cox_normal_lib* result, best_box_cox_normal_lib* bootstrap_results) {

	hist_builder<double> hist(decimals);
	for (int i = 0; i < n; i++)
		hist.add(input[i]);

	best_box_cox_normal<double> internal_result;
	best_box_cox_normal<double>* internal_bootstrap_results;
	if (bootstrap > 0)
		internal_bootstrap_results = new best_box_cox_normal<double>[bootstrap];
	else
		internal_bootstrap_results = NULL;

	algorithm_settings settings = { t1min, t1max, t2min, t2max, sd, tol };
	int return_code = 0;
	std::string cost_func(_cost_func);
	if(cost_func == "classic")
		return_code = ri_estimator<double, kosmic::cost_evaluator_ks_classic<double>>::run(hist, settings, internal_result, internal_bootstrap_results, bootstrap, bootstrap_seed, threads);
	else if(cost_func == "ks_log")
		return_code = ri_estimator<double, kosmic::cost_evaluator_ks_log<double>>::run(hist, settings, internal_result, internal_bootstrap_results, bootstrap, bootstrap_seed, threads);
	else
		return_code = 3;

	// Copy results:
	if(return_code == 0) {
		*result = {
			internal_result.lambda,
			internal_result.mu,
			internal_result.sigma,
			internal_result.interval.cost,
			hist.x(internal_result.interval.t1_i),
			hist.x(internal_result.interval.t2_i)
		};
		for(int i = 0; i < bootstrap; i++)
			bootstrap_results[i] = {
				internal_bootstrap_results[i].lambda,
				internal_bootstrap_results[i].mu,
				internal_bootstrap_results[i].sigma,
				internal_bootstrap_results[i].interval.cost,
				hist.x(internal_bootstrap_results[i].interval.t1_i),
				hist.x(internal_bootstrap_results[i].interval.t2_i)
		};
	}

	if (internal_bootstrap_results)
		delete[] internal_bootstrap_results;

	return return_code;
}
