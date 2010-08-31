#ifndef XLINKSCORER_H_
#define XLINKSCORER_H_
#include "objects.h"
#include "MatchCandidate.h"

class XLinkScorer {
 protected:
  SPECTRUM_T* spectrum_;
  SCORER_T* scorer_;
  ION_CONSTRAINT_T* ion_constraint_;
  MatchCandidate* candidate_;
  int charge_;
  ION_SERIES_T* ion_series_;
  
 public:
  XLinkScorer();
  XLinkScorer(SPECTRUM_T* spectrum, int charge);
  virtual ~XLinkScorer();

  FLOAT_T scoreCandidate(MatchCandidate* candidate);


};


#endif
