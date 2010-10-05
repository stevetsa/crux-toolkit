#include "MPSM_Match.h"
#include "MPSM_MatchCollection.h"
#include "MPSM_Scorer.h"

#include <iostream>

#include "DelimitedFile.h"

using namespace std;

bool MPSM_MatchCompare(MATCH_T* m1, MATCH_T* m2) {
  int c1 = get_match_charge(m1);
  int c2 = get_match_charge(m2);

  if (c1 == c2) {

    PEPTIDE_T* p1 = get_match_peptide(m1);
    char* s1 = get_peptide_unshuffled_modified_sequence(p1);
    string string_s1(s1);
    free(s1);
    PEPTIDE_T* p2 = get_match_peptide(m2);
    char* s2 = get_peptide_unshuffled_modified_sequence(p2);
    string string_s2(s2);
    free(s2);
    return string_s1 < string_s2;
  } else {

    return c1 < c2;
  }
}


SCORER_TYPE_T sort_mode;

MPSM_Match::MPSM_Match() {
  parent_ = NULL;
  invalidate();
}

MPSM_Match::MPSM_Match(MATCH_T* match) {
  parent_ = NULL;
  addMatch(match);
  
  rtime_max_diff_ = 0;

}

/*
MPSM_Match::MPSM_Match(const MPSM_Match& mpsm_match) {
  cout <<"Inside MPSM_Match constructor"<<endl;
  parent_ = mpsm_match.parent_;
  for (int idx = 0;idx < mpsm_match.matches_.size(); idx++) {
    addMatch(mpsm_match.matches_[idx]);
  }
}
*/

MPSM_Match::~MPSM_Match() {
  matches_.clear();
}

void MPSM_Match::setParent(MPSM_MatchCollection* parent) {
  this -> parent_ = parent;
}


MPSM_MatchCollection* MPSM_Match::getParent() {
  return parent_;
}

BOOLEAN_T MPSM_Match::addMatch(MATCH_T* match) {

  //check to make sure match doesn't already exist here.
  for (int idx=0;idx<matches_.size();idx++) {
    if (matches_[idx] == match)
      return FALSE;
  }
  matches_.push_back(match);
  sort(matches_.begin(), matches_.end(), MPSM_MatchCompare);

  


  invalidate();
  return TRUE;
}

BOOLEAN_T MPSM_Match::isDecoy() {
  //If one of the matches is null, then this is a null mpsm.
  for (int idx = 0; idx < matches_.size(); idx++) {
    if (get_match_null_peptide(matches_[idx])) {
      return TRUE;
    }
  }
  return FALSE;
}


BOOLEAN_T MPSM_Match::hasRTime(int match_idx) {
  return FALSE;//return has_rtime[match_idx];
}

FLOAT_T MPSM_Match::getRTime(int match_idx) {
  return 0;//return rtimes[match_idx];
}

void MPSM_Match::setRTime(int match_idx, FLOAT_T rtime) {
  //rtimes[match_idx] = rtime;
}

void MPSM_Match::invalidate() {

  xcorr_score_valid_ = false;
  charge_valid_ = false;
  /*
  for (int idx=0;idx < numMatches();idx++) {
    if (rtimes.size() <= idx) {
      rtimes.push_back(0);
      has_rtime.push_back(false);
    } else {
      has_rtime[idx] = false;
    }
  }
  */
}


ChargeIndex& MPSM_Match::getChargeIndex() {
  if (!charge_valid_) {
    //TODO: validate the charge_index.
    charge_index_.clear();
    for (int idx=0;idx < numMatches();idx++) {
      int charge = get_match_charge(getMatch(idx));
      charge_index_.add(charge);
    }
    charge_valid_ = true;
  }

  return charge_index_;
}

bool MPSM_Match::isChargeHomogeneous() {
  return getChargeIndex().isHomogeneous();
}


MATCH_T* MPSM_Match::getMatch(int match_idx) {
  return matches_[match_idx];
}

MATCH_T* MPSM_Match::operator[] (int match_idx) {
  return matches_[match_idx];
}

int MPSM_Match::numMatches() {
  return matches_.size();
}


FLOAT_T MPSM_Match::getScore(SCORER_TYPE_T match_mode) {

  if (!xcorr_score_valid_) {
    FLOAT_T score = MPSM_Scorer::scoreMPSM(*this, XCORR);
    setScore(match_mode, score);
    return score;
  }
  return xcorr_score_;
}

void MPSM_Match::setScore(SCORER_TYPE_T match_mode, FLOAT_T score) {

  if (match_mode != XCORR) {
    carp(CARP_FATAL,"Score type not implemented for mpsms!");
  }

  xcorr_score_ = score;
  xcorr_score_valid_ = true;

}

void MPSM_Match::getSpectrumNeutralMasses(vector<FLOAT_T>& neutral_masses) {

  neutral_masses.clear();
  MATCH_T* first_match = getMatch(0);
  SPECTRUM_T* spectrum = get_match_spectrum(first_match);

  ChargeIndex& charge_index = getChargeIndex();

  for (int i=0;i<charge_index.size();i++) {
    FLOAT_T mass = get_spectrum_neutral_mass(spectrum, charge_index[i]);
    neutral_masses.push_back(mass);
  }
}

void MPSM_Match::getPeptideMasses(vector<FLOAT_T>& peptide_masses) {
  peptide_masses.clear();
  MATCH_T* match = NULL;

  for (int idx=0;idx<numMatches();idx++) {
    MATCH_T* match = getMatch(idx);
    PEPTIDE_T* peptide = get_match_peptide(match);
    FLOAT_T peptide_mass = get_peptide_peptide_mass(peptide);
    peptide_masses.push_back(peptide_mass);
  }
}

void MPSM_Match::getPeptideModifiedSequences(vector<string>& modified_sequences) {
  modified_sequences.clear();

  MATCH_T* match = NULL;

  for (int idx = 0; idx < numMatches(); idx++) {
    MATCH_T* match = getMatch(idx);
    PEPTIDE_T* peptide = get_match_peptide(match);
    char* seq = get_peptide_unshuffled_modified_sequence(peptide);
    string string_seq(seq);
    modified_sequences.push_back(string_seq);
    free(seq);
  }
}

void MPSM_Match::setRTimeMaxDiff(FLOAT_T rtime_max_diff) {
  rtime_max_diff_ = rtime_max_diff;
}

FLOAT_T MPSM_Match::getRTimeMaxDiff() {
  return rtime_max_diff_;
}

string MPSM_Match::getString() {
  string ans;

  MATCH_T* match = NULL;
  for (int idx = 0; idx < numMatches(); idx++) {
    MATCH_T* match = getMatch(idx);
    PEPTIDE_T* peptide = get_match_peptide(match);
    char* seq = get_peptide_modified_sequence_with_symbols(peptide);
    string string_seq(seq);
    if (idx == 0) 
      ans += string_seq;
    else
      ans += "," + string_seq;
    free(seq);
  }
  return ans;

}

bool MPSM_Match::operator < (const MPSM_Match& match_obj) const {

  MPSM_Match match1 = *this;
  MPSM_Match match2 = match_obj;

  if (match1.numMatches() != match2.numMatches()) {
    return match1.numMatches() < match2.numMatches();
  }

  if (match1.getChargeIndex() != match2.getChargeIndex()) {
    return match1.getChargeIndex() < match2.getChargeIndex();
  }

  vector<string> mod_seq1;
  match1.getPeptideModifiedSequences(mod_seq1);

  vector<string> mod_seq2;
  match2.getPeptideModifiedSequences(mod_seq2);

  for (int idx=0;idx < match1.numMatches();idx++) {
      if (mod_seq1[idx] != mod_seq2[idx]) {
        return (mod_seq1[idx] < mod_seq2[idx]);
      } 
  }
  return false;

}

bool MPSM_Match::operator ==(MPSM_Match& match_obj) {
  if (numMatches() != match_obj.numMatches()) return false;
  if (getChargeIndex() != match_obj.getChargeIndex()) return false;

  vector<string> mod_seq1;
  getPeptideModifiedSequences(mod_seq1);

  vector<string> mod_seq2;
  match_obj.getPeptideModifiedSequences(mod_seq2);

  for (int idx = 0; idx < match_obj.numMatches(); idx++) {
    if (mod_seq1[idx] != mod_seq2[idx]) {
      return false;
    }
  }
  return true;
}


double MPSM_Match::getDeltaCN() {
  return (parent_ -> calcDeltaCNMatch(getScore(XCORR)));
}

double MPSM_Match::getZScore() {
  double xcorr = getScore(XCORR);
  return (parent_ -> calcZScore(xcorr));
}

double MPSM_Match::getSpectrumRTime() {
  return (parent_ -> getSpectrumRTime());
}

double MPSM_Match::getPredictedRTime() {
  return (parent_ -> getPredictedRTime(*this));
}

void MPSM_Match::setXCorrRank(int xcorr_rank) {
  xcorr_rank_ = xcorr_rank;
}

int MPSM_Match::getXCorrRank() {
  if (numMatches() == 1) {
    return get_match_rank(getMatch(0), XCORR);
  } else {
    return xcorr_rank_;
  }
}

string MPSM_Match::getSRankString() {

  ostringstream oss;

  oss << get_match_rank(getMatch(0), XCORR);

  for (int idx=1;idx<numMatches();idx++) {
    oss << "," << get_match_rank(getMatch(idx), XCORR);
  }

  return oss.str();

}


ostream& operator <<(ostream& os, MPSM_Match& match_obj) {

  //cout <<"operator start"<<endl;
  int precision = get_int_parameter("precision");
  os << setprecision(precision);

  //cout <<"getting first_match"<<endl;
  MATCH_T* first_match = match_obj.getMatch(0);
  //cout <<"Getting parent"<<endl;
  MPSM_MatchCollection* parent = match_obj.getParent();
  if (parent == NULL) {
    carp(CARP_WARNING, "Match parent is null!FIX");
  }

  //cout <<"Getting spectrum"<<endl;
  SPECTRUM_T* spectrum = get_match_spectrum(first_match);
  //cout <<"getting scan"<<endl;
  int scan = get_spectrum_first_scan(spectrum);
  
  vector<FLOAT_T> spectrum_neutral_masses;
  match_obj.getSpectrumNeutralMasses(spectrum_neutral_masses);
  
  vector<FLOAT_T> peptide_masses;
  match_obj.getPeptideMasses(peptide_masses);

  ChargeIndex& charge_index = match_obj.getChargeIndex();

  FLOAT_T spectrum_precursor_mz = get_spectrum_precursor_mz(spectrum);

  FLOAT_T xcorr_score = match_obj.getScore(XCORR);
  //int xcorr_rank = match_obj.getScore(XCORR_RANK);
  
  int matches_spectrum = 0;
  if (parent != NULL) {
    matches_spectrum = parent -> numMatches();
  }
  //cout<<"Nummatches:"<<matches_spectrum<<endl;

  vector<string> peptide_modified_sequences;
  match_obj.getPeptideModifiedSequences(peptide_modified_sequences);
  


  

  os << scan << "\t"
     << charge_index << "\t"
     << spectrum_precursor_mz << "\t"
     << DelimitedFile::splice(spectrum_neutral_masses, ',') << "\t" //spectrum neutral mass
     << DelimitedFile::splice(peptide_masses, ',') << "\t" //peptide mass
     << match_obj.getDeltaCN() << "\t" //delta_cn
     /*<< "TODO"*/ << "\t" //sp score
     /*<< "TODO"*/ << "\t" //sp rank
     << xcorr_score << "\t" //xcorr score
     << match_obj.getXCorrRank() << "\t" //xcorr rank.
     /*<< "TODO"*/ << "\t" //p-value
     << match_obj.getZScore() << "\t" //Weibull est. q-value
     /*<< "TODO"*/ << "\t" //decoy q-value (xcorr)
     /*<< "TODO"*/ << "\t" //percolator score
     /*<< "TODO"*/ << "\t" //percolator rank
     /*<< "TODO"*/ << "\t" //percolator q-value
     /*<< "TODO"*/ << "\t" //q-ranker score
     /*<< "TODO"*/ << "\t" //q-ranker q-value
     /*<< "TODO"*/ << "\t" //b/y ions matched
     /*<< "TODO"*/ << "\t" //b/y ions total
     << matches_spectrum << "\t" //matches_spectrum
     << DelimitedFile::splice(peptide_modified_sequences, ',') << "\t" //sequence
     /*<< "TODO"*/ << "\t" //cleavage type
     /*<< "TODO"*/ << "\t" //protein_id
     /*<< "TODO"*/ << "\t" //flanking aa
     /*<< "TODO"*/ << "\t" //unshuffled sequence
     << match_obj.getSpectrumRTime() << "\t" //eta (spectrum retention time)
     << match_obj.getPredictedRTime() << "\t" //beta (predicted retention time)
     << match_obj.getSRankString() << "\t" //shift (spsm ranks)
     << match_obj.getRTimeMaxDiff() /*<< "TODO"*/; //corr (rtime max diff)
     ;

  //cout <<"operator done"<<endl;
  return os;
}

bool compareMPSM_Match(const MPSM_Match& c1, const MPSM_Match& c2) {
  return c1.xcorr_score_ > c2.xcorr_score_;
}


void MPSM_Match::sortMatches(vector<MPSM_Match>& matches, SCORER_TYPE_T match_mode) {
  sort_mode = match_mode;
  sort(matches.begin(), matches.end(), compareMPSM_Match);
}

