#ifndef KROKHINRETENTIONPREDICTOR_H
#define KROKHINRETENTIONPREDICTOR_H

//As implemented in Krokhin et al. 2004.

#include "RetentionPredictor.h"

#include "Match.h"

#include <map>

class KrokhinRetentionPredictor: public RetentionPredictor {
  protected:
    std::map<char, double> aa_rc_coef;
    std::map<char, double> aa_rc_nt_coef;
    double slope;
    double intercept;
  public:
    KrokhinRetentionPredictor();

    virtual ~KrokhinRetentionPredictor();

    virtual FLOAT_T predictRTime(Match* match);

    FLOAT_T predictRTimeS(const char* sequence);
    FLOAT_T predictRTimeS(const char* sequence, int N);
};



#endif
