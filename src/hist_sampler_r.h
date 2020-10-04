#pragma once

#include "kosmic_input.h"

#include <Rcpp.h>
#include <algorithm>
#include <vector>

namespace kosmic {

/*
 Creates random samples from a given histogram using Rcpp::sample()
 */
template<typename T> class hist_sampler_r {
  
  int class_count;
  int entries;
  Rcpp::NumericVector* probabilities;
  
public:
  
  hist_sampler_r(const hist_builder<T>& hist) {
    class_count = hist.classes();
    entries = hist.n();
    probabilities = new Rcpp::NumericVector(class_count);
    for (int i = 0; i < class_count; i++)
      (*probabilities)[i] = (double) hist[i] / entries;
  }
  
  ~hist_sampler_r() {
    delete probabilities;
  }
  
  /*
   Create cdf.
   */
  void cdf(int seed, int temp_counts[], double cdf[]) {
    std::fill(&temp_counts[0], &temp_counts[class_count], 0);
    
    // Sample
    Rcpp::IntegerVector s = Rcpp::sample(class_count, entries, true,
                             (*probabilities), false);
    
    // Count frequency in sample
    for (int i = 0; i < entries; i++)
      temp_counts[s[i]]++;
    
    // Generate a CDF
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