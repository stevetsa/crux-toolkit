#include "AKlammerStaticRetentionPredictor.h"
#include "AKlammerWeights.h"

AKlammerStaticRetentionPredictor::AKlammerStaticRetentionPredictor() {

  weights_composition=(double*)&yeast_20cm_weights_composition;
  weights_cterm=(double*)&yeast_20cm_weights_cterm;
  weights_nterm=(double*)&yeast_20cm_weights_nterm;
  weights_KorR=yeast_20cm_weights_KorR;
  weights_mass=yeast_20cm_weights_mass;
  weights_length=yeast_20cm_weights_length;

  slope = 1;
  intercept = 0;

}


AKlammerStaticRetentionPredictor::~AKlammerStaticRetentionPredictor() {
}


FLOAT_T AKlammerStaticRetentionPredictor::predictRTime(MATCH_T* match) {

  PEPTIDE_T* peptide = get_match_peptide(match);
  int length = get_peptide_length(peptide);

  char* sequence = get_peptide_sequence_pointer(peptide);

  FLOAT_T mass = get_peptide_peptide_mass(peptide);

  double ans = 0.0;

  //COMPOSITION FEATURES
  for (int seq_idx = 0; seq_idx < length; seq_idx++) {
  
    int aa_idx = map_aa_to_idx[sequence[seq_idx]-'A'];
    if (aa_idx == -1) {
      carp(CARP_FATAL,"Cannot find weight for %c",sequence[seq_idx]);
    }

    ans += weights_composition[aa_idx];
  }

  //NTERMINUS FEATURES
  int nterm_idx = map_aa_to_idx[sequence[0]-'A'];

  if (nterm_idx == -1) {
    carp(CARP_FATAL,"Cannot find nterm weight for %c",sequence[0]);
  }

  ans += weights_nterm[nterm_idx];

  //CTERMINUS FEATURES
  int cterm_idx = map_aa_to_idx[sequence[length-2]-'A'];

  if (cterm_idx == -1) {
    carp(CARP_FATAL,"Cannot find cterm weight for %c",sequence[length-2]);
  }

  ans += weights_cterm[cterm_idx];

  //KorR
  if (sequence[length-1] == 'R') {
    ans += weights_KorR;
  } else if (sequence[length-1] != 'K') {
    carp(CARP_DEBUG, "cterm is not K or R! %c", sequence[length-1]);
    //carp(CARP_DEBUG, "%s",sequence);

    //carp(CARP_FATAL,"Exiting...");
  }

  //LENGTH
  ans+= (double)length * weights_length;
  
  //MASS
  ans+= mass * weights_mass;

  //carp(CARP_INFO,"Ans:%lf",ans);
  //carp(CARP_INFO,"rtime:%lf",(slope*ans + intercept));
  return slope * ans + intercept;
  
}

