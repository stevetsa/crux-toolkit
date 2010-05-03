#include "RetentionPredictor.h"
#include "PalmbaldRetentionPredictor.h"
#include "KrokhinRetentionPredictor.h"

#include "parameter.h"


using namespace std;

RetentionPredictor::RetentionPredictor() {
}

RetentionPredictor::~RetentionPredictor() {
}


FLOAT_T RetentionPredictor::predictRTime(MATCH_T* match) {
  /*override this function*/
  return 0.0;
}

double RetentionPredictor::calcMaxDiff(MPSM_Match& mpsm_match) {
  if (mpsm_match.numMatches() <= 1) return 0.0;

  for (int match_idx=0;
    match_idx < mpsm_match.numMatches();
    match_idx++) {
      if (!mpsm_match.hasRTime(match_idx)) {
        mpsm_match.setRTime(match_idx, 
          predictRTime(mpsm_match[match_idx]));
      }
  }

  double max_diff = 0;
  
  for (int idx1 = 0;idx1 < mpsm_match.numMatches()-1;idx1++) {
    for (int idx2 = idx1 + 1;idx2 < mpsm_match.numMatches();idx2++) {

      double current_diff = (mpsm_match.getRTime(idx1) -
        mpsm_match.getRTime(idx2));

      if (fabs(current_diff) > fabs(max_diff)) {
        max_diff = current_diff;
      }
    }
  }

  return max_diff;

    
}


RetentionPredictor* RetentionPredictor::createRetentionPredictor() {

  RTP_TYPE_T rtp_type = get_rtp_type_parameter("rtime-predictor");

  switch(rtp_type) {
    case RTP_KROKHIN:
      carp(CARP_DEBUG,"creating krokhin retention predictor");
      return new KrokhinRetentionPredictor();
    case RTP_PALMBALD:
      carp(CARP_DEBUG,"creating palmbald retention predictor");
      return new PalmbaldRetentionPredictor();
    case RTP_INVALID:
    default:
      carp(CARP_FATAL,"Invalid retention time predictor");
      
  }
}
