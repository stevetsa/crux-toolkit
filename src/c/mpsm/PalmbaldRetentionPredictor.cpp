#include "PalmbaldRetentionPredictor.h"

using namespace std;


PalmbaldRetentionPredictor::PalmbaldRetentionPredictor() {

  aa_coef.insert(make_pair('A', 0.41)); //Alanine
  aa_coef.insert(make_pair('R',-0.76)); //Arganine
  aa_coef.insert(make_pair('N',-0.54)); //Asparagine
  aa_coef.insert(make_pair('D', 0.04)); //Aspartic Acid
  aa_coef.insert(make_pair('C', 1.32)); //Cystiene
  aa_coef.insert(make_pair('E',-0.26)); //Glutamic acid
  aa_coef.insert(make_pair('Q', 1.02)); //Glutamine
  aa_coef.insert(make_pair('G', 0.29)); //Glycine
  aa_coef.insert(make_pair('H', 0.57)); //Histidine
  aa_coef.insert(make_pair('I', 2.70)); //Isoleucine
  aa_coef.insert(make_pair('L', 2.28)); //Leucine
  aa_coef.insert(make_pair('K',-0.66)); //Lysine
  aa_coef.insert(make_pair('M', 0.98)); //Methionine
  aa_coef.insert(make_pair('F', 2.68)); //Phenylalanine
  aa_coef.insert(make_pair('P', 0.97)); //Proline
  aa_coef.insert(make_pair('S',-0.71)); //Serine
  aa_coef.insert(make_pair('T', 0.37)); //Threonine
  aa_coef.insert(make_pair('W', 4.68)); //Tryptophan
  aa_coef.insert(make_pair('Y', 2.78)); //Tyrosine
  aa_coef.insert(make_pair('V', 2.44)); //Valine

  t0 = 4.77;

}

PalmbaldRetentionPredictor::~PalmbaldRetentionPredictor() {
}

FLOAT_T PalmbaldRetentionPredictor::predictRTime(MATCH_T* match) {

  PEPTIDE_T* peptide = get_match_peptide(match);
  int sequence_length = get_peptide_length(peptide);
  char* sequence = get_peptide_sequence_pointer(peptide);

  double ans = t0;

  for (int idx=0;idx<sequence_length;idx++) {

    ans += aa_coef[sequence[idx]];
  }

  return ans;
}
