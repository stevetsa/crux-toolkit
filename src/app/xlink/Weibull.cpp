#include "Weibull.h"
#include "util/crux-utils.h"
#include "util/GlobalParams.h"
#include "model/Scorer.h"
static const FLOAT_T MIN_XCORR_SHIFT = -5.0;
static const FLOAT_T MAX_XCORR_SHIFT  = 5.0;
//#define CORR_THRESHOLD 0.995   // Must achieve this correlation, else punt.
static const FLOAT_T CORR_THRESHOLD = 0.5;
static const FLOAT_T XCORR_SHIFT = 0.05;

using namespace std;

Weibull::Weibull() {
    reset();
}

Weibull::~Weibull() {
    //
}

void Weibull::reset() {
  correlation_ = 0;
  sequences_.clear();
  scores_.clear();
  fit_called_ = false;
  fit_success_ = false;
  duplicates_ = 0;
}

FLOAT_T Weibull::getEta() {
    return eta_;
}

FLOAT_T Weibull::getBeta() {
    return beta_;
}

FLOAT_T Weibull::getShift() {
    return shift_;
}

FLOAT_T Weibull::getCorrelation() {
    return correlation_;
}


void Weibull::addPoint(
  const std::string& sequence,
  FLOAT_T score
) {

  set<string>::iterator it = sequences_.find(sequence);

  if (it == sequences_.end()) {
    scores_.push_back(score);
    sequences_.insert(sequence);
    fit_called_ = false;
  } else {
    duplicates_++;
  }

}

void Weibull::fit() {
    
  int nscores = scores_.size();
  
  carp(CARP_DEBUG, "Fitting %i scores", nscores);
  carp(CARP_DEBUG, "There were %i dups detected during insert", duplicates_);
  
    
  sort(scores_.begin(), scores_.end(), greater<FLOAT_T>());
  int num_tail_samples = (int)(nscores * GlobalParams::getFractionToFit());
  carp(CARP_DEBUG, "num tail:%d", num_tail_samples);  
  FLOAT_T max_shift = scores_[0];
  FLOAT_T min_shift = min(MIN_XCORR_SHIFT, scores_[nscores-1]);
    
  if (nscores > 29) {
    FLOAT_T* scores_ptr = scores_.data();
    max_shift = scores_ptr[29];
    fit_three_parameter_weibull(
      scores_ptr,
      num_tail_samples,
      nscores,
      min_shift,
      max_shift,
      XCORR_SHIFT,
      0,
      &eta_,
      &beta_,
      &shift_,
      &correlation_
    );
    if (correlation_ < CORR_THRESHOLD) {
      carp(CARP_WARNING, "Weibull fit failed corr:%g n:%d nt:%d "
        "min:%g max:%g mins:%g maxs:%g",
        correlation_,
        nscores,
        num_tail_samples,
        scores_[nscores-1],
        scores_[0],
        min_shift,
        max_shift
      );
      fit_success_ = false;
    } else {
      fit_success_ = true;
    }  
  } else {
    carp(CARP_DEBUG, "Too few scores to fit weibull");
    fit_success_ = false;    
  }
  //Fit was called.
  fit_called_ = true;
}

FLOAT_T Weibull::getWeibullPValue(FLOAT_T score, bool logp) {
  if (!fit_called_) {
    fit();
  }
  
  FLOAT_T pvalue = compute_weibull_pvalue(score, eta_, beta_, shift_);
  if (logp) {
    pvalue = log10(pvalue);
  }
  return(pvalue);
  
}

FLOAT_T Weibull::getECDFPValue(FLOAT_T score, bool logp) {
  if (!fit_called_) {
    fit();
  }
  
  FLOAT_T pvalue;
  
  //end point check, this is probably unstable but a conservative estimate.
  if (score > scores_[0]) {
    pvalue = 1.0 / (1.0 + (FLOAT_T)scores_.size());
  } else {
    //TODO, use upper/lower bound?
    int idx = 0;
    while(idx < scores_.size() && score <= scores_[idx]) {
        idx++;
    }
    pvalue = (FLOAT_T)idx / (FLOAT_T)scores_.size();
  }
  
  if (logp) {
    pvalue = log10(pvalue);
  }
  return(pvalue);
  
  
}

FLOAT_T Weibull::getPValue(FLOAT_T score, bool logp) {
  if (!fit_called_) {
    fit();
  }
  FLOAT_T pvalue;
  if (fit_success_) {
    pvalue = getWeibullPValue(score, logp);
  }
  //If the fit failed or the estimated p-value is invalid, then
  //return the ECDF
  if (!fit_success_ || pvalue == 0 || pvalue != pvalue) {
    //x != x should test for NaN.
    if (pvalue == 0 || pvalue != pvalue) {
      carp(CARP_WARNING, "invalid pvalue! %g, returning ecdf", pvalue);
    }
    pvalue = getECDFPValue(score, logp);
  }
  
  return(pvalue);
}




