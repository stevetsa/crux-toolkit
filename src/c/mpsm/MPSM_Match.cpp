#include "MPSM_Match.h"
#include "MPSM_MatchCollection.h"
#include "MPSM_Scorer.h"

#include "SpectrumZState.h"

#include <iostream>

#include "DelimitedFile.h"

using namespace std;


int compareSequences(MATCH_T* m1, MATCH_T* m2) {
  PEPTIDE_T* p1 = get_match_peptide(m1);
  PEPTIDE_T* p2 = get_match_peptide(m2);

  int l1 = get_peptide_length(p1);
  int l2 = get_peptide_length(p2);

  char* s1 = get_peptide_sequence_pointer(p1);
  char* s2 = get_peptide_sequence_pointer(p2);

  int idx = 0;

  while (idx < l1 && idx < l2) {
  
    if (s1[idx] < s2[idx]) {
      return -1;
    } else if (s1[idx] > s2[idx]) {
      return 1;
    }
    //otherwise they are equal
    idx++;
  }
  
  if (l1 == l2) {
    return 0;
  } else if (l1 < l2) {
    return -1;
  } else {
    return 1;
  }


}

bool MPSM_MatchCompare(MATCH_T* m1, MATCH_T* m2) {


  if (get_match_zstate(m1) != get_match_zstate(m2)) {
    return get_match_zstate(m1) < get_match_zstate(m2);
  } else {
    return (compareSequences(m1, m2) == -1);
  }
}



bool MATCH_T_Compare2(MATCH_T* m1, MATCH_T* m2) {
  if (get_match_charge(m1) != get_match_charge(m2)) {
    return get_match_charge(m1) < get_match_charge(m2);
  } else {
    return (compareSequences(m1, m2) == -1);
  }
}

bool compareMPSM_MatchVisited(const MPSM_Match& m1, const MPSM_Match& m2) {

  if (m1.numMatches() != m2.numMatches()) {
    return m1.numMatches() < m2.numMatches();
  }

  vector<MATCH_T*> matches_1(m1.matches_.begin(), m1.matches_.end());
  
  vector<MATCH_T*> matches_2(m2.matches_.begin(), m2.matches_.end());

  sort(matches_1.begin(), matches_1.end(), MATCH_T_Compare2);      
  sort(matches_2.begin(), matches_2.end(), MATCH_T_Compare2);

  for (int idx=0;idx<matches_1.size();idx++) {
      MATCH_T* submatch_m1 = matches_1.at(idx);
      MATCH_T* submatch_m2 = matches_2.at(idx);

      int charge1 = get_match_charge(submatch_m1);
      int charge2 = get_match_charge(submatch_m2);

      if (charge1 != charge2) {
        return (charge1 < charge2);
        break;
      }

      int compare = compareSequences(submatch_m1, submatch_m2);

      if (compare != 0) {
        return (compare == -1); //return whether the sequence in m1 < m2.
      }

  }

  return false; //all charges/sequences are equal, so not <.
  
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

bool isMATCH_T_Equal(MATCH_T* m1, MATCH_T* m2) {
  if (get_match_charge(m1) != get_match_charge(m2)) {
    return false;
  }
  PEPTIDE_T* p1 = get_match_peptide(m1);
  char* s1 = get_peptide_unshuffled_modified_sequence(p1);
  string string_s1(s1);
  free(s1);
  
  PEPTIDE_T* p2 = get_match_peptide(m2);
  char* s2 = get_peptide_unshuffled_modified_sequence(p2);
  string string_s2(s2);
  free(s2);

  return string_s1 == string_s2;
}


BOOLEAN_T MPSM_Match::addMatch(MATCH_T* match) {

  //check to make sure match doesn't already exist here.
  for (int idx=0;idx<matches_.size();idx++) {
    if (isMATCH_T_Equal(matches_[idx], match)) {
      return false;
    }  
  }
  matches_.push_back(match);
  sort(matches_.begin(), matches_.end(), MPSM_MatchCompare);

  


  invalidate();
  return true;
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

  zstate_valid_ = false;

  for (int idx=0;idx < NUMBER_SCORER_TYPES;idx++) {
    have_match_score_[idx] = false;
    have_match_rank_[idx] = false;
  }
}


ZStateIndex& MPSM_Match::getZStateIndex() {
  if (!zstate_valid_) {
    //TODO: validate the charge_index.
    zstate_index_.clear();
    for (int idx=0;idx < numMatches();idx++) {
      MATCH_T* match = getMatch(idx);
      zstate_index_.add(get_match_zstate(match));
    }
    zstate_valid_ = true;
  }

  return zstate_index_;
}
/*
bool MPSM_Match::isChargeHomogeneous() {
  return getChargeIndex().isHomogeneous();
}
*/




MATCH_T* MPSM_Match::getMatch(int match_idx) const {
  return matches_.at(match_idx);
}

MATCH_T* MPSM_Match::operator[] (int match_idx) {
  return matches_[match_idx];
}

int MPSM_Match::numMatches() const {
  return matches_.size();
}

Spectrum* MPSM_Match::getSpectrum() {

  return get_match_spectrum(getMatch(0));
}

char* MPSM_Match::getSequence(int match_idx) {
  MATCH_T* match = getMatch(match_idx);
  PEPTIDE_T* peptide = get_match_peptide(match);

  return get_peptide_sequence(peptide);
}

MODIFIED_AA_T* MPSM_Match::getModifiedSequence(int match_idx) {
  return get_peptide_modified_aa_sequence(get_match_peptide(getMatch(match_idx)));
} 

int MPSM_Match::getMaxCharge() {

  int max_charge = getCharge(0);

  for (int idx=1;idx < numMatches();idx++) {
    max_charge = max(max_charge, getCharge(idx));
  }

  return max_charge;

}

int MPSM_Match::getCharge(int match_idx) {
  get_match_charge(getMatch(match_idx));
}


FLOAT_T MPSM_Match::getScore(SCORER_TYPE_T match_mode) {

  if (!have_match_score_[match_mode]) {

    setScore(match_mode, MPSM_Scorer::score(*this, XCORR));
  }

  return match_scores_[match_mode];
}

FLOAT_T MPSM_Match::getScoreConst(SCORER_TYPE_T match_mode) const {

  if (!have_match_score_[match_mode]) {
    carp(CARP_FATAL," Dont have score for %d",(int)match_mode);
  }

  return match_scores_[match_mode];

}

void MPSM_Match::setScore(SCORER_TYPE_T match_mode, FLOAT_T score) {

  match_scores_[match_mode] = score;
  have_match_score_[match_mode] = true;

}


void MPSM_Match::setRank(SCORER_TYPE_T match_mode, int rank) {
  
  match_rank_[match_mode] = rank; 
  have_match_rank_[match_mode] = true;
}

int MPSM_Match::getRank(SCORER_TYPE_T match_mode) {
  if (!have_match_rank_[match_mode]) {

    carp(CARP_FATAL,"Don't have match rank for %d!",(int)match_mode);

  }

  return match_rank_[match_mode];
}


void MPSM_Match::getSpectrumNeutralMasses(vector<FLOAT_T>& neutral_masses) {

  neutral_masses.clear();
  MATCH_T* first_match = getMatch(0);
  Spectrum* spectrum = get_match_spectrum(first_match);

  ZStateIndex& zstate_index = getZStateIndex();

  for (int i=0;i<zstate_index.size();i++) {
    FLOAT_T mass = zstate_index[i].getNeutralMass();
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

FLOAT_T MPSM_Match::getRTimeAveDiff() {
  
  return 0;

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

  if (match1.getZStateIndex() != match2.getZStateIndex()) {
    return match1.getZStateIndex() < match2.getZStateIndex();
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
  if (getZStateIndex() != match_obj.getZStateIndex()) return false;

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

string MPSM_Match::getSRankString() {

  ostringstream oss;

  oss << get_match_rank(getMatch(0), XCORR);

  for (int idx=1;idx<numMatches();idx++) {
    oss << "," << get_match_rank(getMatch(idx), XCORR);
  }

  return oss.str();

}

int MPSM_Match::getFirstScan() {

  return getSpectrum()->getFirstScan();
}

string MPSM_Match::getChargeString() {
  return getZStateIndex().getChargeString();
}

FLOAT_T MPSM_Match::getSpectrumPrecursorMZ() {
  return getSpectrum()->getPrecursorMz();

}

string MPSM_Match::getPeptideMassString() {

  vector<FLOAT_T> peptide_masses;
  getPeptideMasses(peptide_masses);

  return DelimitedFile::splice(peptide_masses, ',');
  
}

string MPSM_Match::getNeutralMassString() {

  vector<FLOAT_T> neutral_masses;
  getSpectrumNeutralMasses(neutral_masses);
  
  return DelimitedFile::splice(neutral_masses, ',');

}

string MPSM_Match::getSequenceString() {
  
  vector<string> sequences;
  getPeptideModifiedSequences(sequences);

  return DelimitedFile::splice(sequences, ',');

}

int MPSM_Match::getMatchesPerSpectrum() {
  
  return getParent() -> numMatches();

}

FLOAT_T MPSM_Match::getXCorrSumDiff() {
  
  FLOAT_T ans = getScore(XCORR);

  for (int idx=0;idx < matches_.size();idx++) {

    ans -= get_match_score(matches_[idx], XCORR);
  }

  return ans;

}

FLOAT_T MPSM_Match::getXCorrMaxDiff() {

  FLOAT_T max_score = get_match_score(matches_[0], XCORR);

  for (int idx=1;idx < matches_.size();idx++) {
    FLOAT_T current_score = get_match_score(matches_[idx], XCORR);
    if (current_score > max_score) {
      max_score = current_score;
    }
  }

  return getScore(XCORR) - max_score;
}

FLOAT_T MPSM_Match::getAreaRatio() {

  if (matches_.size() == 1) { 
    return 1.0;
  }

  double sum = 0.0;

  double min_area = -1;

  bool min_unset = true;

  for (unsigned int idx=0;idx<matches_.size();idx++) {
    
    SpectrumZState& current_zstate = 
      get_match_zstate(matches_[idx]);

    FLOAT_T area = current_zstate.getArea();

    if ((min_unset) || (area < min_area)) {
      min_area = area;
      min_unset = false;
    }

    sum += area;
  }

  if (sum == 0.0) return 1;

  return min_area / sum * (FLOAT_T)matches_.size();
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
  Spectrum* spectrum = get_match_spectrum(first_match);
  //cout <<"getting scan"<<endl;
  int scan = spectrum->getFirstScan();
  
  vector<FLOAT_T> spectrum_neutral_masses;
  match_obj.getSpectrumNeutralMasses(spectrum_neutral_masses);
  
  vector<FLOAT_T> peptide_masses;
  match_obj.getPeptideMasses(peptide_masses);

  ZStateIndex& zstate_index = match_obj.getZStateIndex();

  FLOAT_T spectrum_precursor_mz = spectrum->getPrecursorMz();

  FLOAT_T xcorr_score = 0;
  if (get_boolean_parameter("mpsm-do-sort")) {
    xcorr_score = match_obj.getScore(XCORR);
  } else {
    xcorr_score = 0;//match_obj.xcorr_score_;
  }
  //int xcorr_rank = match_obj.getScore(XCORR_RANK);
  
  int matches_spectrum = 0;
  if (parent != NULL) {
    matches_spectrum = parent -> numMatches();
  }
  //cout<<"Nummatches:"<<matches_spectrum<<endl;

  vector<string> peptide_modified_sequences;
  match_obj.getPeptideModifiedSequences(peptide_modified_sequences);
  


  

  os << scan << "\t"
     << zstate_index << "\t"
     << spectrum_precursor_mz << "\t"
     << DelimitedFile::splice(spectrum_neutral_masses, ',') << "\t" //spectrum neutral mass
     << DelimitedFile::splice(peptide_masses, ',') << "\t" //peptide mass
     << match_obj.getDeltaCN() << "\t" //delta_cn
     /*<< "TODO"*/ << "\t" //sp score
     /*<< "TODO"*/ << "\t" //sp rank
     << xcorr_score << "\t" //xcorr score
     << match_obj.getRank(XCORR) << "\t" //xcorr rank.
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

  return c1.getScoreConst(sort_mode) > c2.getScoreConst(sort_mode);
}


void MPSM_Match::sortMatches(vector<MPSM_Match>& matches, SCORER_TYPE_T match_mode) {
  sort_mode = match_mode;
  sort(matches.begin(), matches.end(), compareMPSM_Match);
}

bool CompareMPSM_MatchVisited::operator() (
  const MPSM_Match& m1, 
  const MPSM_Match& m2) const {

    return compareMPSM_MatchVisited(m1, m2);
}
