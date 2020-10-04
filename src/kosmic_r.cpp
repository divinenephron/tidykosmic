
#include "kosmic.h"
#include "hist_sampler_r.h"

#include <Rcpp.h>
using namespace Rcpp;

//' Call Kosmic algorithm
//' 
// [[Rcpp::export]]
List kosmic_impl(NumericVector input_vector, int decimals, int bootstrap, int bootstrap_seed, int threads,
                              double t1min, double t1max, double t2min, double t2max, double sd, double tol){
  int n = input_vector.size();
  
  kosmic::hist_builder<double> hist(decimals);
  for (int i = 0; i < n; i++){
    hist.add(input_vector[i]);
  }
  
  kosmic::best_box_cox_normal<double> internal_result;
  kosmic::best_box_cox_normal<double>* internal_bootstrap_results;
  
  if(bootstrap > 0){
    internal_bootstrap_results = new kosmic::best_box_cox_normal<double>[bootstrap];
  }else{
    internal_bootstrap_results = NULL;
  }
  kosmic::algorithm_settings settings = { t1min, t1max, t2min, t2max, sd, tol };
  
  int return_code = kosmic::ri_estimator<double, kosmic::cost_evaluator_ks_classic<double>>::run(hist, settings, internal_result, internal_bootstrap_results, bootstrap, bootstrap_seed, threads);
  // Copy results:
  if(return_code != 0) {
    return List::create(Named("return_code", return_code));
  }
  
  NumericVector r_result = NumericVector::create(
    internal_result.lambda,
    internal_result.mu,
    internal_result.sigma,
    internal_result.interval.cost,
    hist.x(internal_result.interval.t1_i),
    hist.x(internal_result.interval.t2_i)
  );
  
  NumericVector r_boot = NumericVector(Dimension(bootstrap, 6));
  for(int i = 0; i < bootstrap; i++) {
    r_boot(i, 0) = internal_bootstrap_results[i].lambda;
    r_boot(i, 1) = internal_bootstrap_results[i].mu;
    r_boot(i, 2) = internal_bootstrap_results[i].sigma;
    r_boot(i, 3) = internal_bootstrap_results[i].interval.cost;
    r_boot(i, 4) = hist.x(internal_bootstrap_results[i].interval.t1_i);
    r_boot(i, 5) = hist.x(internal_bootstrap_results[i].interval.t2_i);
  }
  
  return List::create(
    Named("return_code", return_code),
    Named("result", r_result),
    Named("boot", r_boot));
}

//' Resample data using the same code as the Kosmic agorithm
//' 
// [[Rcpp::export]]
NumericVector kosmic_resample(NumericVector results, NumericVector counts, int decimals) {
  // Make histogram
  kosmic::hist_builder<double> hist(decimals);
  for (int i = 0; i < results.size(); i++) {
    for (int j = 0; j < counts[i]; j++) {
      hist.add(results[i]);
    }
  }
  
  // Resample
  int* temp_counts = new int[hist.classes()];
  double* temp_cdfs = new double[hist.classes()];
  kosmic::hist_sampler_r<double> hist_sampler(hist);
  hist_sampler.cdf(0, temp_counts, temp_cdfs);
  
  // Return resampled frequencies
  NumericVector res(hist.classes());
  for (int i = 0; i < hist.classes(); i++) {
    res[i] = temp_counts[i];
  }
  return res;
}
