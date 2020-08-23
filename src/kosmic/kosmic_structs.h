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

namespace kosmic {

	/*
	 kosmic algorithm settings.
	*/
	struct algorithm_settings {
		double t1min_p, t1max_p;
		double t2min_p, t2max_p;
		double sd_p;
		double abs_tol;
	};

	/*
	 Bounds and penalty of an evaluated truncation interval, x >= t1, x < t2.
	*/
	template<typename T> struct evaluated_interval {
		T cost;
		int t1_i, t2_i;

		evaluated_interval() {}
		evaluated_interval(T cost) : cost(cost) {}
		evaluated_interval(T cost, int t1_i, int t2_i) : cost(cost), t1_i(t1_i), t2_i(t2_i) {}

		T value() const { return cost; }

	};

	template<typename T> struct best_normal {
		T mu, sigma;
		evaluated_interval<T> interval;
	};

	template<typename T> struct best_box_cox_normal {
		T lambda, mu, sigma;
		evaluated_interval<T> interval;
	};

}

