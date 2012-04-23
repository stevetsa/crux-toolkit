#include "XLinkMatch.h"

#include "parameter.h"
#include "Scorer.h"
#include <sstream>
#include <ios>
#include <iomanip>
#include <iostream>


using namespace std;

XLinkMatch::XLinkMatch() {
  parent_ = NULL;
  pvalue_= 1;
  for (int idx = 0;idx < NUMBER_MASS_TYPES;idx++) {
    mass_calculated_[idx] = false;
    mass_[idx] = 0;
  }
}

XLinkMatch::~XLinkMatch() {

}

void XLinkMatch::decrementPointerCount() {

  pointer_count_--;
}

void XLinkMatch::computeWeibullPvalue(
  FLOAT_T shift,
  FLOAT_T eta,
  FLOAT_T beta) {

  pvalue_ = compute_weibull_pvalue(getScore(XCORR), eta, beta, shift);
  cerr <<"eta:"<<eta<<" beta:"<<beta<<" shift:"<<shift<<endl;
  cerr <<"xcorr:"<<getScore(XCORR)<<" pvalue:"<<pvalue_<<endl;
}
/*
void XLinkMatch::setBYIonsMatched(int by_ions_matched) {
  by_ions_matched_ = by_ions_matched;
}

int XLinkMatch::getBYIonsMatched() {
  return by_ions_matched_;
}

void XLinkMatch::setBYIonsTotal(int by_ions_total) {
  by_ions_total_ = by_ions_total;
}

int XLinkMatch::getBYIonsTotal() {
  return by_ions_total_;
}
*/
string XLinkMatch::getProteinIdString() {
  Peptide* peptide = this -> getPeptide(0);

  if (peptide == NULL) {
    return string("");
  } else {
    return XLink::get_protein_ids_locations(peptide);
  }
}


FLOAT_T XLinkMatch::getMass(MASS_TYPE_T mass_type) {

  if (!mass_calculated_[mass_type]) {
    mass_[mass_type] = calcMass(mass_type);
    mass_calculated_[mass_type] = true;
  }
  return mass_[mass_type];
}


FLOAT_T XLinkMatch::getPPMError() {
  FLOAT_T mono_mass = getMass(MONO);
  
  return (mono_mass - parent_->getSpectrumNeutralMass()) / mono_mass * 1e6;
  
}

void XLinkMatch::setParent(XLinkMatchCollection* parent) {
  parent_ = parent;
}

/**
 * Print one field in the tab-delimited output file, based on column index.
 */
void XLinkMatch::printOneMatchField(
  int      column_idx,             ///< Index of the column to print. -in
  MatchCollection* collection,  ///< collection holding this match -in 
  MatchFileWriter*    output_file,            ///< output stream -out
  int      scan_num,               ///< starting scan number -in
  FLOAT_T  spectrum_precursor_mz,  ///< m/z of spectrum precursor -in
  int      num_target_matches,            ///< target matches for this spectrum -in
  int      num_decoy_matches, ///< decoy matches for this spectrum -in
  int      b_y_total,              ///< total b/y ions -in
  int      b_y_matched             ///< Number of b/y ions matched. -in
) {

  //cerr <<"XLinkMatch::printOneMatchField:";
  //cerr <<get_column_header(column_idx)<<endl;
  switch ((MATCH_COLUMNS_T)column_idx) {

  case PEPTIDE_MASS_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, getMass(MONO));
    break;
  case PVALUE_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, pvalue_);
    break;
  case SEQUENCE_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
      getSequenceString());
    break;
  case PROTEIN_ID_COL:
    output_file->setColumnCurrentRow(
      (MATCH_COLUMNS_T)column_idx, 
      getProteinIdString()); 
    break;
  case FLANKING_AA_COL:
    break;
  default:
    Match::printOneMatchField(column_idx,
      collection,
      output_file,
      scan_num,
      spectrum_precursor_mz,
      num_target_matches,
      num_decoy_matches,
      b_y_total,
      b_y_matched
    );
  }
}
