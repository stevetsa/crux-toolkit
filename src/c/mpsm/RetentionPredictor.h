#ifndef RETENTIONPREDICTOR_H_
#define RETENTIONPREDICTOR_H_

#include "Match.h"
#include "MPSM_Match.h"


class RetentionPredictor {

  public:
    RetentionPredictor();
    
    virtual ~RetentionPredictor();

    //override this.
    virtual FLOAT_T predictRTime(Match* match);
    
    //create an instance of the retention predictor, depending
    //upon the parameter.
    static RetentionPredictor* createRetentionPredictor();

    double calcMaxDiff(MPSM_Match& mpsm_match);

};


#endif
