#ifndef MPSM_SCORER_H
#define MPSM_SCORER_H
#include "objects.h"

#include "MPSM_Match.h"


class MPSM_Scorer {
  protected:
    SPECTRUM_T* spectrum_;
    int max_charge_;
    SCORER_T* scorer_;
    FLOAT_T* theoretical_;
    int max_mz_;

    ION_CONSTRAINT_T** ion_constraints_;

  public:

    MPSM_Scorer(SPECTRUM_T* spectrum, int max_charge);
    virtual ~MPSM_Scorer();

    FLOAT_T calcScore(MPSM_Match& mpsm_match, SCORER_TYPE_T match_mode);

    static FLOAT_T scoreMPSM(MPSM_Match& mpsm_match, SCORER_TYPE_T match_mode); 
};

#endif
