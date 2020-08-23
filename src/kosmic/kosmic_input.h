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
#include <vector>
#include <random>

namespace kosmic {

	/*
	 Histogram builder class.
	 Does not store actual values but their counts using an internal std::vector.
	*/
	template<typename T> class hist_builder {

		std::vector<int> counts;
		int entries;
		T first, last, factor;

	public:

		hist_builder() : entries(0), first(INFINITY), last(INFINITY), factor(INFINITY) {}

		hist_builder(int decimals) : entries(0), first(INFINITY), last(INFINITY) {
			set_decimals(decimals);
		}

		void set_decimals(int decimals) {
			assert(entries == 0);
			assert(decimals >= 0 && decimals <= 9);
			factor = std::pow((T) 10.0, decimals);
		}

		/*
		 Add value to histogram.
		*/
		void add(T value) {
			assert(value >= 0.0);
			int index = (int) std::round(value * factor - first);
			if (index >= 0 && (unsigned) index < counts.size())
				counts[index]++;
			else {
				value = std::round(value * factor);
				if (counts.size() == 0) {
					first = value;
					last = value;
					counts.push_back(1);
				} else if (value < first) {
					int n = (int) std::round(first - value);
					counts.resize(n + counts.size());
					std::move(&counts[0], &counts[counts.size() - n], &counts[n]);
					counts[0] = 1;
					for (int i = 1; i < n; i++)
						counts[i] = 0;
					first = value;
				} else if (value > last) {
					int n = (int) std::round(value - last);
					counts.resize(n + counts.size());
					for (int i = (int) counts.size() - n;  i < (int) counts.size() - 1; i++)
						counts[i] = 0;
					counts[counts.size() - 1] = 1;
					last = value;
				}
			}
			entries++;
		}

		/*
		 Count of classes in histogram.
		*/
		int classes() const {
			return (int) counts.size();
		}

		/*
		 Number of values in histogram.
		*/
		int n() const {
			return entries;
		}

		/*
		 First/minimum class.
		*/
		T min() const {
			return first / factor;
		}

		/*
		 Last/maximum class.
		*/
		T max() const {
			return last / factor;
		}

		/*
		 Count of values with given index in histogram.
		*/
		int operator[](int index) const {
			return counts[index];
		}

		/*
		 Return x position at given index.
		*/
		T x(int index) const {
			if(counts.size() > 1)
				return (first + (T)index * (last - first) / (counts.size() - (int)1)) / factor;
			else
				return first / factor;
		}

	};

	/*
	 Creates random (Mersenne twister/std::mt19937) samples (sample with replacement) from a given histogram.
	*/
	template<typename T> class hist_sampler {

		int class_count;
		int entries;
		double* probabilities;
		std::discrete_distribution<int>* dist;

	public:

		hist_sampler(const hist_builder<T>& hist) {
			class_count = hist.classes();
			entries = hist.n();
			probabilities = new double[class_count];
			for (int i = 0; i < class_count; i++)
				probabilities[i] = (double) hist[i] / entries;
			dist = new std::discrete_distribution<int>(probabilities, &probabilities[class_count]);
		}

		~hist_sampler() {
			delete[] probabilities;
			delete dist;
		}

		/*
		 Create cdf.
		*/
		void cdf(int seed, int temp_counts[], double cdf[]) {
			std::fill(&temp_counts[0], &temp_counts[class_count], 0);
			std::mt19937 mt(seed);
			for (int i = 0; i < entries; i++)
				temp_counts[(*dist)(mt)]++;
			int total = 0;
			for (int i = 0; i < class_count; i++) {
				total += temp_counts[i];
				cdf[i] = (double) total / entries;
			}
		}

		/*
		 Count of classes in histogram.
		*/
		int classes() const {
			return class_count;
		}

		/*
		 Number of values in histogram.
		*/
		int n() const {
			return entries;
		}

	};

	/*
	 Cumulative density "function": array of calculated cumulative densities.
	*/
	template<typename T> struct cdf {
	protected:

		void set_array(T values[], int length, int decimals) {
			// Round & sort values:
			T factor = std::pow((T) 10.0, decimals);
			for (int i = 0; i < length; i++)
				values[i] = std::round(values[i] * factor) / factor;
			std::sort(values, &values[length]);

			// Calculate x positions:
			T delta = std::pow((T) 10.0, -decimals);
			classes = (int) (((values[length - 1]) - values[0]) / delta + 1);
			x = new T[classes];
			for (int i = 0; i < classes; i++)
				x[i] = values[0] + ((T) i / (classes - 1)) * (values[length - 1] - values[0]);

			// Calculate cumulative frequency:
			cum_freq = new T[classes];
			int gr_i = 0;
			for (int i = 0; i < classes; i++) {
				while (values[gr_i] <= x[i] && gr_i < length)
					gr_i++;
				cum_freq[i] = (T) gr_i / length;
			}
		}
		/*
		void set_const_array(const T values[], int length, int decimals) {
			T* values_temp = new T[length];
			std::copy(values, &values[length], values_temp);
			set_array(values_temp, length, decimals);
			delete[] values_temp;
		}*/

		void set_hist(const hist_builder<T>& hist) {
			classes = hist.classes();
			x = new T[classes];
			for (int i = 0; i < classes; i++)
				x[i] = hist.min() + ((T) i * ((hist.max() - hist.min()) / (classes - (int) 1)));
			cum_freq = new T[classes];
			int total = 0;
			for (int i = 0; i < classes; i++) {
				total += hist[i];
				cum_freq[i] = (T) total / hist.n();
			}
		}

		void free() {
			if (x) {
				delete[] x, cum_freq;
				x = cum_freq = NULL;
			}
		}

	public:
		int classes;
		T* x;
		T* cum_freq;

		cdf() : classes(0), x(NULL), cum_freq(NULL) {}
		cdf(const hist_builder<T>& hist) { set_hist(hist); }
		~cdf() { free(); }

	};

};
