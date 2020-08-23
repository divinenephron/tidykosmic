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

#include <cassert>
#include <cmath>

#include "../lib/erfinv.h"

#include "kosmic_structs.h"

namespace kosmic {

	/*
	 CDF (% of values <= x) of a normal distribution defined by mu/sigma.
	*/
	template<typename T> T std_normal_cdf(T mu, T sigma, T x) {
		const T SQRT_2 = 1.4142135623730950488016887242097;
		return 0.5* (1.0 + erf((x - mu) / (sigma * SQRT_2)));
	}

	/*
	 Quantile of a normal distribution defined by mu/sigma.
	*/
	template<typename T> T quantile_std_normal(T mu, T sigma, T q) {
		const T SQRT_2 = 1.4142135623730950488016887242097;
		return mu + sigma * SQRT_2 * erfinv(2 * q - 1);
	}

	/*
	 Box-Cox transformation of x.
	*/
	template<typename T> T box_cox_transform(T x, T lambda) {
		if (lambda == 0.0)
			return std::log(x);
		else
			return (std::pow(x, lambda) - 1) / lambda;
	}

	/*
	 Inverse Box-Cox transformation of x.
	*/
	template<typename T> T box_cox_inverse(T x, T lambda) {
		if (lambda == 0.0)
			return std::exp(x);
		else
			return pow(x * lambda + 1.0, 1.0 / lambda);
	}

	/*
	 Returns the quantile q (0.0-1.0) of an array of sorted values.
	*/
	template<typename T> T quantile(T q, const T sorted_values[], int length) {
		assert(q >= 0 && q <= 1);
		assert(length > 0);
		T index = q * (length-1);
		T index_integer = std::floor(index);
		T index_frac = index - index_integer;
		if (index_frac == 0)
			return sorted_values[(int) index_integer];
		else
			return sorted_values[(int) index_integer] + index_frac * (sorted_values[(int) index_integer + 1] - sorted_values[(int) index_integer]);
	}

	/*
	 Returns the quantile q (0.0-1.0) of a Box-Cox transformed distribution.
	*/
	template<typename T> T quantile(T q, const best_box_cox_normal<T>& params) {
		return box_cox_inverse(quantile_std_normal<T>(params.mu, params.sigma, q), params.lambda);
	}

}
