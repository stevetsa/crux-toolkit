#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <limits>
#include <vector>

int carp_threshold = 0;

#define CARP_DETAILED_DEBUG 2
#define CARP_DEBUG 1

void carp(int level, const char* format, ...) {
  if (carp_threshold < level)
    return;
  va_list ap;
  va_start(ap, format);
  vprintf(format, ap);
  va_end(ap);
}

/**
 * Fits a two-parameter Weibull distribution to the input data. 
 *
 * Called by the three parameter weibull fitting function to see if
 * the proposed shift gives the best correlation.  If there are too
 * few data points, sets correlation to 0 (minimum value).
 * http:// www.chinarel.com/onlincebook/LifeDataWeb/rank_regression_on_y.htm
 * \returns eta, beta and the correlation coefficient.
 */
void fit_two_parameter_weibull(
    double* data, ///< the data to be fit. should be in descending order -in
    int fit_data_points, ///< the number of data points to fit -in
    int total_data_points, ///< the total number of data points -in
    double shift, ///< the amount by which to shift our data -in
    double* eta,      ///< the eta parameter of the Weibull dist -out
    double* beta,      ///< the beta parameter of the Weibull dist -out
    double* correlation ///< the best correlation -out
    ){

  double* X = (double*)malloc(sizeof(double) * fit_data_points); //hold data here

  // transform data into an array of values for fitting
  // shift (including only non-neg data values) and take log
  int idx;
  for(idx=0; idx < fit_data_points; idx++){
    double score = data[idx] + shift; // move right by shift
    if (score <= 0.0){
      carp(CARP_DEBUG, "Reached negative score at idx %i", idx);
      fit_data_points = idx;
      break;
    } 
    X[idx] = log(score);
    // carp(CARP_DEBUG, "X[%i]=%.6f=ln(%.6f)", idx, X[idx], score);
  }

  double* Y   = (double*)malloc(sizeof(double) * fit_data_points);
  for(idx=0; idx < fit_data_points; idx++){
    int reverse_idx = total_data_points - idx;
    double F_T_idx = (reverse_idx - 0.3) / (total_data_points + 0.4);
    Y[idx] = log( -log(1.0 - F_T_idx) );
    //carp(CARP_DEBUG, "Y[%i]=%.6f", idx, Y[idx]);
  }

  int N = fit_data_points; // rename for formula's sake
  double sum_Y  = 0.0;
  double sum_X  = 0.0;
  double sum_XY = 0.0;
  double sum_XX = 0.0;
  for(idx=0; idx < fit_data_points; idx++){
    sum_Y  += Y[idx];
    sum_X  += X[idx];
    sum_XX += X[idx] * X[idx];
    sum_XY += X[idx] * Y[idx];
  }
  carp(CARP_DETAILED_DEBUG, "sum_X=%.6f", sum_X);
  carp(CARP_DETAILED_DEBUG, "sum_Y=%.6f", sum_Y);
  carp(CARP_DETAILED_DEBUG, "sum_XX=%.6f", sum_XX);
  carp(CARP_DETAILED_DEBUG, "sum_XY=%.6f", sum_XY);

  double b_num    = sum_XY - (sum_X * sum_Y / N);
  carp(CARP_DETAILED_DEBUG, "b_num=%.6f", b_num);
  double b_denom  = sum_XX - sum_X * sum_X / N;
  carp(CARP_DETAILED_DEBUG, "b_denom=%.6f", b_denom);
  double b_hat    = b_num / b_denom;

  double a_hat    = (sum_Y - b_hat * sum_X) / N;
  *beta = b_hat;
  *eta  = exp( - a_hat / *beta );

  double c_num   = 0.0;
  double c_denom_X = 0.0;
  double c_denom_Y = 0.0;
  double mean_X = sum_X / N;
  double mean_Y = sum_Y / N;
  for (idx=0; idx < N; idx++){
    double X_delta = X[idx] - mean_X; 
    double Y_delta = Y[idx] - mean_Y;
    c_num += X_delta * Y_delta;
    c_denom_X += X_delta * X_delta;
    c_denom_Y += Y_delta * Y_delta;
  }
  double c_denom = sqrt(c_denom_X * c_denom_Y);
  if (c_denom == 0.0){
    //carp(CARP_FATAL, "Zero denominator in correlation calculation!");
    carp(CARP_DETAILED_DEBUG, "Zero denominator in correlation calculation!");
    *correlation = 0.0; // min value
    *eta = 0;
    *beta = 0;
  } else {
    *correlation = c_num / c_denom;
  }
  carp(CARP_DETAILED_DEBUG, "shift=%.6f",shift);
  carp(CARP_DETAILED_DEBUG, "eta=%.6f", *eta);
  carp(CARP_DETAILED_DEBUG, "beta=%.6f", *beta);
  carp(CARP_DETAILED_DEBUG, "correlation=%.6f", *correlation);

  free(Y);
  free(X);
}

/**
 * Fits a three-parameter Weibull distribution to the input data. 
 * Implementation of Weibull distribution parameter estimation from 
 * http:// www.chinarel.com/onlincebook/LifeDataWeb/rank_regression_on_y.htm
 * \returns eta, beta, c (which in this case is the amount the data should
 * be shifted by) and the best correlation coefficient
 */
void fit_three_parameter_weibull(
    double* data, ///< the data to be fit -in
    int fit_data_points, ///< the number of data points to fit -in
    int total_data_points, ///< the total number of data points to fit -in
    double min_shift, ///< the minimum shift to allow -in
    double max_shift, ///< the maximum shift to allow -in
    double step,      ///< step for shift -in
    double corr_threshold, ///< minimum correlation, else no fit -in
    double* eta,      ///< the eta parameter of the Weibull dist -out
    double* beta,      ///< the beta parameter of the Weibull dist -out
    double* shift,     ///< the best shift -out
    double* correlation   ///< the best correlation -out
    ){
  
  double correlation_tolerance = 0.1;
  
  double best_eta = 0.0;
  double best_beta = 0.0;
  double best_shift = 0.0;
  double best_correlation = 0.0;

  double cur_eta = 0.0;
  double cur_beta = 0.0;
  double cur_correlation = 0.0;
  double cur_shift = 0.0;

  for (cur_shift = max_shift; cur_shift > min_shift ; cur_shift -= step){

    fit_two_parameter_weibull(data, fit_data_points, total_data_points, 
			      cur_shift, &cur_eta, &cur_beta, &cur_correlation);

    if (cur_correlation > best_correlation){
      best_eta = cur_eta;
      best_beta = cur_beta;
      best_shift = cur_shift;
      best_correlation = cur_correlation;
    } else if (cur_correlation < best_correlation - correlation_tolerance){
      break;
    }
  }

  // Only store the parameters if the fit was good enough.
  *correlation = best_correlation;
  if (best_correlation >= corr_threshold) {
    *eta = best_eta;
    *beta = best_beta;
    *shift = best_shift;
  } else {
    *eta = 0.0;
    *beta = 0.0;
    *shift = 0.0;
  }
}


/**
 * Compute a p-value for a given score w.r.t. a Weibull with given parameters.
 *\returns the p_value
 */
double compute_weibull_pvalue(
  double score, ///< The score for the scoring peptide -in
  double eta,   ///< The eta parameter of the Weibull -in
  double beta,  ///< The beta parameter of the Weibull -in
  double shift  ///< The shift parameter of the Weibull -in
  ){
  carp(CARP_DETAILED_DEBUG, "Stat: score = %.6f", score);

  double return_value;

  // No Weibull parameter, return NaN.
  if (eta == 0.0) {
    carp(CARP_DETAILED_DEBUG, "Failed fit, returning p-value=NaN");
    return_value = std::numeric_limits<double>::quiet_NaN();
  }
  // undefined past shift, give lowest possible score.
  else if (score + shift <= 0) {
    carp(CARP_DETAILED_DEBUG, "Bad shift, returning p-value=1");
    return_value = 1.0;
  }
  else {
    return_value = exp(-pow((score + shift) / eta, beta));
    carp(CARP_DETAILED_DEBUG, "Stat: pvalue before = %g", return_value);
  }
  return(return_value);
}

#define MIN_WEIBULL_MATCHES 40
#define MIN_SHIFT -5.0
#define MAX_SHIFT  5.0
//#define CORR_THRESHOLD 0.995   // Must achieve this correlation, else punt.
#define CORR_THRESHOLD 0.0       // For now, turn off the threshold.
#define SHIFT_STEP 0.05

double GetPValue(double* scores, int num_scores, int num_candidates) {
  double eta, beta, shift, correlation;
  if (num_candidates < MIN_WEIBULL_MATCHES)
    return std::numeric_limits<double>::quiet_NaN();
  if (num_scores > num_candidates * 0.55)
    num_scores = int(num_candidates * 0.55);
  fit_three_parameter_weibull(scores+1, num_scores-1, num_candidates,
			      MIN_SHIFT, MAX_SHIFT, SHIFT_STEP, CORR_THRESHOLD,
			      &eta, &beta, &shift, &correlation);
  return compute_weibull_pvalue(scores[0], eta, beta, shift);
}
