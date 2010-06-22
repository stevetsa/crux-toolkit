#ifndef KROKHINRETENTIONPREDICTOR_H
#define KROKHINRETENTIONPREDICTOR_H

#include "RetentionPredictor.h"

#include "match.h"

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

    virtual FLOAT_T predictRTime(MATCH_T* match);
    FLOAT_T predictRTimeS(const char* sequence);
};



#endif
