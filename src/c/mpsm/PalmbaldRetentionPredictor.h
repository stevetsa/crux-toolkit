#ifndef PALMBALDRETENTIONPREDICTOR_H
#define PALMBALDRETENTIONPREDICTOR_H

#include "RetentionPredictor.h"

#include "match.h"

#include <map>

class PalmbaldRetentionPredictor: public RetentionPredictor {
  protected:
    double t0;
    std::map<char, double> aa_coef;

  public:
    PalmbaldRetentionPredictor();

    virtual ~PalmbaldRetentionPredictor();

    virtual FLOAT_T predictRTime(MATCH_T* match);

};



#endif
