#ifndef NULLRETENTIONPREDICTOR_H
#define NULLRETENTIONPREDICTOR_H

#include "RetentionPredictor.h"

#include "Match.h"

#include <map>

class NullRetentionPredictor: public RetentionPredictor {
  protected:
    std::map<char, double> aa_rc_coef;
    std::map<char, double> aa_rc_nt_coef;
    double slope;
    double intercept;
  public:
    NullRetentionPredictor();

    virtual ~NullRetentionPredictor();

    virtual FLOAT_T predictRTime(Match* match);

    FLOAT_T predictRTimeS(const char* sequence);
    FLOAT_T predictRTimeS(const char* sequence, int N);
};



#endif
