#include "KrokhinRetentionPredictor.h"

using namespace std;


KrokhinRetentionPredictor::KrokhinRetentionPredictor() {

  aa_rc_coef.insert(make_pair('W', 11.0)); //Tryptophan
  aa_rc_coef.insert(make_pair('F', 10.5)); //Phenylalanine
  aa_rc_coef.insert(make_pair('L', 9.6));  //Leucine
  aa_rc_coef.insert(make_pair('I', 8.4));  //Isoleucine
  aa_rc_coef.insert(make_pair('M', 5.8));  //Methionine
  aa_rc_coef.insert(make_pair('V', 5.0));  //Valine
  aa_rc_coef.insert(make_pair('Y', 4.0));  //Tyrosine
  aa_rc_coef.insert(make_pair('A', 0.8));  //Alanine
  aa_rc_coef.insert(make_pair('T', 0.4));  //Threonine
  aa_rc_coef.insert(make_pair('P', 0.2));  //Proline
  aa_rc_coef.insert(make_pair('Q', 0.0));  //Glutamine
  aa_rc_coef.insert(make_pair('D', -0.5)); //Aspartic Acid
  aa_rc_coef.insert(make_pair('C',-0.8));  //Cystiene (carboxamidomethylated)
  aa_rc_coef.insert(make_pair('S',-0.8));  //Serine
  aa_rc_coef.insert(make_pair('E',-0.9));  //Glutamic acid
  aa_rc_coef.insert(make_pair('G',-0.9));  //Glycine  
  aa_rc_coef.insert(make_pair('N',-1.2));  //Asparagine
  aa_rc_coef.insert(make_pair('R',-1.3));  //Arganine
  aa_rc_coef.insert(make_pair('H',-1.3));  //Histidine
  aa_rc_coef.insert(make_pair('K',-1.9)); //Lysine
    
  aa_rc_nt_coef.insert(make_pair('W', -4.0)); //Tryptophan
  aa_rc_nt_coef.insert(make_pair('F', -7.0)); //Phenylalanine
  aa_rc_nt_coef.insert(make_pair('L', -9.0));  //Leucine
  aa_rc_nt_coef.insert(make_pair('I', -8.0));  //Isoleucine
  aa_rc_nt_coef.insert(make_pair('M', -5.5));  //Methionine
  aa_rc_nt_coef.insert(make_pair('V', -5.5));  //Valine
  aa_rc_nt_coef.insert(make_pair('Y', -3.0));  //Tyrosine
  aa_rc_nt_coef.insert(make_pair('A', -1.5));  //Alanine
  aa_rc_nt_coef.insert(make_pair('T', 5.0));  //Threonine
  aa_rc_nt_coef.insert(make_pair('P', 4.0));  //Proline
  aa_rc_nt_coef.insert(make_pair('Q', 7.0));  //Glutamine
  aa_rc_nt_coef.insert(make_pair('D', 9.0)); //Aspartic Acid
  aa_rc_nt_coef.insert(make_pair('C', 4.0));  //Cystiene (carboxamidomethylated)
  aa_rc_nt_coef.insert(make_pair('S', 5.0));  //Serine
  aa_rc_nt_coef.insert(make_pair('E', 1.0));  //Glutamic acid
  aa_rc_nt_coef.insert(make_pair('G', 5.0));  //Glycine  
  aa_rc_nt_coef.insert(make_pair('N', 5.0));  //Asparagine
  aa_rc_nt_coef.insert(make_pair('R', 8.0));  //Arganine
  aa_rc_nt_coef.insert(make_pair('H', 4.0));  //Histidine
  aa_rc_nt_coef.insert(make_pair('K', 4.6)); //Lysine

  slope = 1.0;
  intercept = 0.0;

}

KrokhinRetentionPredictor::~KrokhinRetentionPredictor() {
}

FLOAT_T KrokhinRetentionPredictor::predictRTime(MATCH_T* match) {
  char* sequence = get_match_sequence(match);
  double ans = predictRTimeS(sequence);
  free(sequence);
  return ans;
}

FLOAT_T KrokhinRetentionPredictor::predictRTimeS(const char* sequence) {

  int N = strlen(sequence);
  
  double K_L;

  if (N < 10) {
    K_L = 1.0 - 0.027 * (double)(10 - N);
  } else if (N > 20) {
    K_L = 1.0 - 0.014 * (double)(N - 20);
  } else {
    K_L = 1.0;
  }

  double sum_Rc = 0.0;
  for (int idx=0;idx<N;idx++) {
    sum_Rc += aa_rc_coef.find(sequence[idx]) -> second;
  }

  double Rc_1_Nt = aa_rc_nt_coef.find(sequence[0]) -> second;
  double Rc_2_Nt = aa_rc_nt_coef.find(sequence[1]) -> second;
  double Rc_3_Nt = aa_rc_nt_coef.find(sequence[2]) -> second;

  double H = K_L * (sum_Rc + 0.42*Rc_1_Nt + 0.22*Rc_2_Nt + 0.05*Rc_3_Nt);
  
  

  if (H >= 38.0) {
    H = H - 0.3 * (H - 38.0);
  }
    
  double ans = intercept + H * slope;

  return ans;
}
