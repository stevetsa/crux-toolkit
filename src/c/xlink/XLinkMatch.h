#ifndef XLINKMATCH_H_
#define XLINKMATCH_H_


class XLinkMatch {
 protected:
  XLinkPeptide xlink_peptide_;
  FLOAT_T xcorr_score_;
  FLOAT_T ppm_error_;
  FLOAT_T pvalue_;
};


#endif
