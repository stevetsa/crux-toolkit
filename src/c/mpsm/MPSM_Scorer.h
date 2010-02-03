#ifndef MPSM_SCORER_H
#define MPSM_SCORER_H
#include "objects.h"

#include "MPSM_Match.h"


class MPSM_Scorer {


  public:
    static FLOAT_T scoreMPSM(MPSM_Match& mpsm_match, SCORER_TYPE_T match_mode); 
};

#endif
