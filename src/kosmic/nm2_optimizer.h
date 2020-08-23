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
#include <utility>

/*
 Nelder Mead optimization implementation for 2-arg functions.
 T: precision, e.g. double or float
 Objective_T: type returned by Evaluator_T::evaluate(x0, x1), must implement T value()
 Evaluator_T: object which is evaluated, must implement Objective_T evaluate(T x0, T x1)
*/
template<typename T, class Objective_T, class Evaluator_T> class nm2_optimizer {
protected:

	/*
	 2d point helper class.
	*/
	struct point {
		T x, y;

		point() {}
		point(T x, T y) : x(x), y(y) {}

		point operator+(point right) { return point(x + right.x, y + right.y); }
		point operator-(point right) { return point(x - right.x, y - right.y); }
		point operator*(T right) { return point(x * right, y * right); }
		point operator/(T right) { return point(x / right, y / right); }

	};

	/*
	 2d point, with objective value at that point, helper class.
	*/
	struct evaluated_point {
		point p;
		Objective_T objective;

		evaluated_point() {}
		evaluated_point(T x, T y, Evaluator_T& evaluator) : p(x, y) {
			objective = evaluator.evaluate(p.x, p.y);
		}
		evaluated_point(point p, Evaluator_T& evaluator) : p(p) {
			objective = evaluator.evaluate(p.x, p.y);
		}

	};

	// Simplex:
	evaluated_point s0, s1, s2;

	// Objective function evaluator:
	Evaluator_T& evaluator;

	/*
	 Constructor, creates initial simplex and evaluates the objective function at the initial simplex' points.
	*/
	nm2_optimizer(T x0, T x1, Evaluator_T& evaluator) : evaluator(evaluator) {
		s0 = evaluated_point(x0, x1, evaluator);
		if (s0.p.x != 0.0)
			s1 = evaluated_point(x0 * (T) 1.05, x1, evaluator);
		else
			s1 = evaluated_point((T) 0.00025, x1, evaluator);
		if (s0.p.y != 0.0)
			s2 = evaluated_point(x0, x1 * (T) 1.05, evaluator);
		else
			s2 = evaluated_point(x0, (T) 0.00025, evaluator);
	}

	/*
	 Constructor, creates initial simplex and evaluates the objective function at the initial simplex' points.
	*/
	nm2_optimizer(T x0_0, T x1_0, T x0_1, T x1_1, T x0_2, T x1_2, Evaluator_T& evaluator) : evaluator(evaluator) {
		s0 = evaluated_point(x0_0, x1_0, evaluator);
		s1 = evaluated_point(x0_1, x1_1, evaluator);
		s2 = evaluated_point(x0_2, x1_2, evaluator);
	}

	/*
	 Nelder-Mead algorithm step.
	*/
	Objective_T step() {
		// Sort by objective:
		if (s0.objective.value() < s1.objective.value()) {
			if (s1.objective.value() > s2.objective.value()) {
				std::swap(s1, s2);
				if (s0.objective.value() > s1.objective.value())
					std::swap(s0, s1);
			}
		} else {
			if (s1.objective.value() < s2.objective.value()) {
				std::swap(s0, s1);
				if (s1.objective.value() > s2.objective.value())
					std::swap(s1, s2);
			} else
				std::swap(s0, s2);
		}
		assert(s0.objective.value() <= s1.objective.value());
		assert(s1.objective.value() <= s2.objective.value());

		// Algorithm step:
		const T alpha = 1.0, gamma = 2.0, beta = 0.5, sigma = 0.5;
		point m = (s0.p + s1.p) / 2.0;
		evaluated_point r(m * (1.0 + alpha) - s2.p * alpha, evaluator);
		if (r.objective.value() < s0.objective.value()) {
			evaluated_point e(m * (1.0 + gamma) - s2.p * gamma, evaluator);
			s2 = e.objective.value() < r.objective.value() ? e : r;
			return s0.objective;
		} else if (r.objective.value() < s1.objective.value()) {
			s2 = r;
			return s0.objective;
		}
		if (r.objective.value() > s2.objective.value())
			r = s2;
		evaluated_point c(m * beta + r.p * (1 - beta), evaluator);
		if (c.objective.value() < s2.objective.value()) {
			s2 = c;
			return s0.objective;
		}
		s1 = evaluated_point(s0.p * sigma + s1.p * (1 - sigma), evaluator);
		s2 = evaluated_point(s0.p * sigma + s2.p * (1 - sigma), evaluator);
		return s0.objective;
	}

public:

	/*
	 Optimization result.
	*/
	struct Result {
		T best_x0, best_x1;
		Objective_T best_objective;
	};

	/*
	 Optimization function.
	*/
	template<typename Tolerance_T> static Result optimize(T x0, T x1, Evaluator_T& evaluator, Tolerance_T abs_tol) {
		nm2_optimizer nm(x0, x1, evaluator);
		T last_rs = nm.step().value();
		for(;;) {
			T rs = nm.step().value();
			if (rs != last_rs) {
				if (last_rs - rs < abs_tol)
					return Result{ nm.s0.p.x, nm.s0.p.y, nm.s0.objective };
				last_rs = rs;
			}
		}
	}

	/*
	 Optimization function.
	*/
	template<typename Tolerance_T> static Result optimize(T x0_0, T x1_0, T x0_1, T x1_1, T x0_2, T x1_2, Evaluator_T& evaluator, Tolerance_T abs_tol) {
		nm2_optimizer nm(x0_0, x1_0, x0_0, x1_1, x0_2, x1_2, evaluator);
		T last_rs = nm.step().value();
		for (;;) {
			T rs = nm.step().value();
			if (rs != last_rs) {
				if (last_rs - rs < abs_tol)
					return Result{ nm.s0.p.x, nm.s0.p.y, nm.s0.objective };
				last_rs = rs;
			}
		}
	}

};
