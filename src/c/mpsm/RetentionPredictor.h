#ifndef RETENTIONPREDICTOR_H_
#define RETENTIONPREDICTOR_H_

#include "Match.h"
#include "MPSM_Match.h"


class RetentionPredictor {

  protected:
    static RetentionPredictor* predictor_;

  public:
    RetentionPredictor();
    
    virtual ~RetentionPredictor();

    //override this.
    virtual FLOAT_T predictRTime(Match* match);
    
    //create the static instance of the retention predictor, depending
    //upon the parameter.
    static void createRetentionPredictor();
    static RetentionPredictor* getStaticRetentionPredictor();

    double calcMaxDiff(const MPSM_Match& mpsm_match);

    
    


};


#endif
