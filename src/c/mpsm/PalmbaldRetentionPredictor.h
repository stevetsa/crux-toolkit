#ifndef PALMBALDRETENTIONPREDICTOR_H
#define PALMBALDRETENTIONPREDICTOR_H

#include "RetentionPredictor.h"

#include "Match.h"

#include <map>

class PalmbaldRetentionPredictor: public RetentionPredictor {
  protected:
    double t0;
    std::map<char, double> aa_coef;

  public:
    PalmbaldRetentionPredictor();

    virtual ~PalmbaldRetentionPredictor();

    virtual FLOAT_T predictRTime(Match* match);

};



#endif
