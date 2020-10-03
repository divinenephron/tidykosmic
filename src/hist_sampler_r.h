#pragma once

#include "kosmic_input.h"

#include <Rcpp.h>
#include <cassert>
#include <vector>
#include <random>

namespace kosmic {

/*
 Creates random (Mersenne twister/std::mt19937) samples (sample with replacement) from a given histogram.
 */
template<typename T> class hist_sampler_r {
  
  int class_count;
  int entries;
  double* probabilities;
  std::discrete_distribution<int>* dist;
  
public:
  
  hist_sampler_r(const hist_builder<T>& hist) {
    class_count = hist.classes();
    entries = hist.n();
    probabilities = new double[class_count];
    for (int i = 0; i < class_count; i++)
      probabilities[i] = (double) hist[i] / entries;
    dist = new std::discrete_distribution<int>(probabilities, &probabilities[class_count]);
  }
  
  ~hist_sampler_r() {
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

};