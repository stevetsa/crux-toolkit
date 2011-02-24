#ifndef XLINKSCORER_H_
#define XLINKSCORER_H_
#include "objects.h"
#include "MatchCandidate.h"

class XLinkScorer {
 protected:
  Spectrum* spectrum_;
  SCORER_T* scorer_;
  IonConstraint* ion_constraint_;
  MatchCandidate* candidate_;
  int charge_;
  IonSeries* ion_series_;
  
 public:
  XLinkScorer();
  XLinkScorer(Spectrum* spectrum, int charge);
  virtual ~XLinkScorer();

  FLOAT_T scoreCandidate(MatchCandidate* candidate);


};


#endif
