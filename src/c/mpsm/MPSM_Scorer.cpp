#include "MPSM_Scorer.h"




FLOAT_T MPSM_Scorer::scoreMPSM(MPSM_Match& mpsm_match, SCORER_TYPE_T match_mode) {
  //only allow XCorr to be scored right now.

  if (match_mode != XCORR) {
    carp(CARP_FATAL,"Score type not implemented for mpsms!");
  }

  //Touch right now.
  mpsm_match.getChargeIndex();


}
