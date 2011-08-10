#include "XLinkMatch.h"

#include "parameter.h"
#include "scorer.h"
#include <sstream>
#include <ios>
#include <iomanip>
#include <iostream>


using namespace std;

XLinkMatch::XLinkMatch() {
  parent_ = NULL;
  pvalue_= 1;
  for (int idx = 0;idx < NUMBER_MASS_TYPES;idx++) {
    mass_calculated_[idx] = FALSE;
    mass_[idx] = 0;
  }
}

XLinkMatch::~XLinkMatch() {

}

void XLinkMatch::computeWeibullPvalue(
  FLOAT_T shift,
  FLOAT_T eta,
  FLOAT_T beta) {

  pvalue_ = compute_weibull_pvalue(getScore(XCORR), eta, beta, shift);
}

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

string XLinkMatch::getProteinIdString(int peptide_idx) {
  PEPTIDE_T* peptide = this -> getPeptide(peptide_idx);

  if (peptide == NULL) {
    return string("");
  } else {
    return XLink::get_protein_ids_locations(peptide);
  }
}


FLOAT_T XLinkMatch::getMass(MASS_TYPE_T mass_type) {

  if (!mass_calculated_[mass_type]) {
    mass_[mass_type] = calcMass(mass_type);
    mass_calculated_[mass_type] = TRUE;
  }
  return mass_[mass_type];
}


FLOAT_T XLinkMatch::getPPMError() {
  FLOAT_T mono_mass = getMass(MONO);
  
  return (mono_mass - parent_->getSpectrumNeutralMass()) / mono_mass * 1e6;
  
}

string XLinkMatch::getResultHeader() {

  ostringstream oss;
  oss << "scan" << "\t"
      << "charge" << "\t"
      << "spectrum precursor m/z" << "\t"
      << "spectrum neutral mass" << "\t"
      << "peptide mass mono" << "\t"
      << "peptide mass average" << "\t"
      << "mass error(ppm)" << "\t";
  if (get_boolean_parameter("compute-sp")) {
    oss << "sp score" << "\t"
        << "sp rank" << "\t"
        << "b/y ions matched" << "\t"
        << "b/y ions total" << "\t";
  }
  oss << "xcorr score" << "\t"
      << "xcorr rank" << "\t"
      << "p-value" << "\t"
      << "matches/spectrum" << "\t"
      << "sequence" << "\t"
      << "protein id(loc) 1"<< "\t"
      << "protein id(loc) 2";

  return oss.str();
}

std::string XLinkMatch::getResultString() {
  ostringstream ss;
  int precision = get_int_parameter("precision");
  ss << std::setprecision(precision);

  ss << parent_->getScan() << "\t" //scan
     << parent_->getCharge() << "\t" //charge
     << parent_->getPrecursorMZ() << "\t" //precursor mz
     << parent_->getSpectrumNeutralMass() << "\t"
     << this->getMass(MONO) << "\t"
     << this->getMass(AVERAGE) << "\t"
     << this->getPPMError() << "\t";

  if (get_boolean_parameter("compute-sp")) {
    ss << this->getScore(SP) << "\t"
       << this->getRank(SP) << "\t"
       << this->getBYIonsMatched() << "\t"
       << this->getBYIonsTotal() << "\t";
  }

  ss << this->getScore(XCORR) << "\t"
     << this->getRank(XCORR) << "\t"
     << pvalue_ << "\t"
//     << parent_->getExperimentSize() << "\t"
     << this->getSequenceString() << "\t"
     << this->getProteinIdString(0) << "\t"   //protein id(loc) 1
     << this->getProteinIdString(1);  //protein id(loc) 2
  string out_string = ss.str();
  return out_string;
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
  int      num_matches,            ///< num matches in spectrum -in
  int      b_y_total,              ///< total b/y ions -in
  int      b_y_matched             ///< Number of b/y ions matched. -in
) {

  cerr <<"XLinkMatch::printOneMatchField:";
  cerr <<get_column_header(column_idx)<<endl;
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
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
      getProteinIdString(0)); //TODO - fix this.
    break;
  case FLANKING_AA_COL:
    break;
  default:
    Match::printOneMatchField(column_idx,
      collection,
      output_file,
      scan_num,
      spectrum_precursor_mz,
      num_matches,
      b_y_total,
      b_y_matched
    );
  }
}
