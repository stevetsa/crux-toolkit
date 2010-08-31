#include "MatchCandidate.h"

#include "parameter.h"
#include "scorer.h"
#include <sstream>
#include <ios>
#include <iomanip>

using namespace std;

void MatchCandidate::computeWeibullPvalue(
  FLOAT_T shift,
  FLOAT_T eta,
  FLOAT_T beta) {

  pvalue_ = compute_weibull_pvalue(xcorr_, eta, beta, shift);
}

void MatchCandidate::setXCorr(FLOAT_T xcorr) {
  xcorr_ = xcorr;
}

FLOAT_T MatchCandidate::getXCorr() {
  return xcorr_;
}

string MatchCandidate::getResultHeader() {

  ostringstream oss;
  oss << "scan" << "\t"
      << "charge" << "\t"
      << "spectrum precursor m/z" << "\t"
      << "spectrum neutral mass" << "\t"
      << "peptide mass mono" << "\t"
      << "peptide mass average" << "\t"
      << "mass error(ppm)" << "\t"
      << "xcorr score" << "\t"
      << "xcorr rank" << "\t"
      << "p-value" << "\t"
      << "matches/spectrum" << "\t"
      << "sequence" << "\t"
      << "protein id"<< "\t"
      << "by total" << "\t"
      << "by observable (0-1200)" << "\t"
      << "by observable bin (0-1200)" << "\t"
      << "by observable (0-max)" << "\t"
      << "by obsrevable bin (0-max)" << "\t"
      << "by observed bin" << "\t"
      << "ion current total" << "\t"
      << "ion current observed" << "\t"
      << "ions observable bin (0-1200)";

  return oss.str();
}

std::string MatchCandidate::getResultString() {
  ostringstream ss;
  int precision = get_int_parameter("precision");
  ss << std::setprecision(precision);

  ss << parent_->getScan() << "\t" //scan
     << parent_->getCharge() << "\t" //charge
     << parent_->getPrecursorMZ() << "\t" //precursor mz
     << parent_->getSpectrumNeutralMass() << "\t"
     << this->getMass() << "\t"
     << this->getMass() << "\t"
     << "" /*this->getPPMError(parent_->getSpectrumNeutralMass())*/ << "\t"
     << xcorr_ << "\t" //xcorr score
     << ""/*xcorr_rank_*/ << "\t"
     << pvalue_ << "\t"
     << parent_->size() << "\t"
     << this->getSequenceString() << "\t"
     << "" << "\t"   //protein id
     << "" << "\t"   // by total
     << "" << "\t"   // by observable (0-1200)
     << "" << "\t"   // by observable bin (0-1200)
     << "" << "\t"   // by observable (0-max)
     << "" << "\t"   // by observable bin (0-max)
     << "" << "\t"   // by observed bin
     << "" << "\t"   // ion current total
     << "" << "\t"   // ion current observed
     << "";   // ions observable bin (0-1200)
  string out_string = ss.str();
  return out_string;
}

void MatchCandidate::setParent(MatchCandidateVector* parent) {
  parent_ = parent;
}
