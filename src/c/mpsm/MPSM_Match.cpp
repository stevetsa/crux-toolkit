#include "MPSM_Match.h"
#include "MPSM_MatchCollection.h"
#include "MPSM_Scorer.h"

#include "RetentionPredictor.h"
#include "SpectrumZState.h"

#include <iostream>

#include "DelimitedFile.h"

using namespace std;


int compareSequences(Match* m1, Match* m2) {
  Peptide* p1 = m1->getPeptide();
  Peptide* p2 = m2->getPeptide();

  int l1 = p1->getLength();
  int l2 = p2->getLength();

  char* s1 = p1->getSequencePointer();
  char* s2 = p2->getSequencePointer();

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

bool MPSM_MatchCompare(Match* m1, Match* m2) {


  if (m1->getZState() != m2->getZState()) {
    return m1->getZState() < m2->getZState();
  } else {
    return (compareSequences(m1, m2) == -1);
  }
}



bool Match_Compare2(Match* m1, Match* m2) {
  if (m1->getCharge() != m2->getCharge()) {
    return m1->getCharge() < m2->getCharge();
  } else {
    return (compareSequences(m1, m2) == -1);
  }
}

bool compareMPSM_MatchVisited(
  const vector<Match*>& m1, 
  const vector<Match*>& m2
  ) {

  if (m1.size() != m2.size()) {
    return m1.size() < m2.size();
  }

  vector<Match*> matches_1(m1.begin(), m1.end());
  
  vector<Match*> matches_2(m2.begin(), m2.end());

  sort(matches_1.begin(), matches_1.end(), Match_Compare2);      
  sort(matches_2.begin(), matches_2.end(), Match_Compare2);

  for (int idx=0;idx<matches_1.size();idx++) {
      Match* submatch_m1 = matches_1.at(idx);
      Match* submatch_m2 = matches_2.at(idx);

      int charge1 = submatch_m1->getCharge();
      int charge2 = submatch_m2->getCharge();

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

MPSM_Match::MPSM_Match(Match* match) {
  parent_ = NULL;
  addMatch(match);

}



MPSM_Match::~MPSM_Match() {
  matches_.clear();
}

void MPSM_Match::setParent(MPSM_MatchCollection* parent) {
  this -> parent_ = parent;
}


MPSM_MatchCollection* MPSM_Match::getParent() {
  return parent_;
}

bool isMatch_Equal(Match* m1, Match* m2) {
  if (m1->getCharge() != m2->getCharge()) {
    return false;
  }
  Peptide* p1 = m1->getPeptide();
  char* s1 = p1->getUnshuffledModifiedSequence();
  string string_s1(s1);
  free(s1);
  
  Peptide* p2 = m2->getPeptide();
  char* s2 = p2->getUnshuffledModifiedSequence();
  string string_s2(s2);
  free(s2);

  return string_s1 == string_s2;
}


BOOLEAN_T MPSM_Match::addMatch(Match* match) {

  //check to make sure match doesn't already exist here.
  for (int idx=0;idx<matches_.size();idx++) {
    if (isMatch_Equal(matches_[idx], match)) {
      return false;
    }  
  }
  matches_.push_back(match);
  sort(matches_.begin(), matches_.end(), MPSM_MatchCompare);

  
  parent_ = NULL;

  invalidate();
  return true;
}

bool MPSM_Match::isDecoy() {
  //If one of the matches is null, then this is a null mpsm.
  for (int idx = 0; idx < matches_.size(); idx++) {
    if (matches_[idx]->getNullPeptide()) {
      return true;
    }
    if (getProteinIdsString().find("random_") != string::npos) {
      return true;
    }
  }
  return false;
}

void MPSM_Match::invalidate() {

  if (parent_ != NULL) {
    parent_->invalidate();
  }
  match_scores_.clear();
  match_rank_.clear();
}


ZStateIndex& MPSM_Match::getZStateIndex() {


  if (parent_ == NULL) {
    carp(CARP_FATAL, "Parent not set!");
  } 

  return parent_ -> getZStateIndex();
}

ZStateIndex MPSM_Match::calcZStateIndex() {
  ZStateIndex zstate_index;
  for (int idx=0;idx < numMatches();idx++) {
    Match* match = getMatch(idx);
    zstate_index.add(match->getZState());
  }

  return zstate_index;
  
}


Match* MPSM_Match::getMatch(int match_idx) const {
  return matches_.at(match_idx);
}

Match* MPSM_Match::operator[] (int match_idx) {
  return matches_[match_idx];
}

int MPSM_Match::numMatches() const {
  return matches_.size();
}

Spectrum* MPSM_Match::getSpectrum() {

  return getMatch(0)->getSpectrum();
}

char* MPSM_Match::getSequence(int match_idx) {
  Match* match = getMatch(match_idx);
  Peptide* peptide = match->getPeptide();

  return peptide->getSequence();
}

MODIFIED_AA_T* MPSM_Match::getModifiedSequence(int match_idx) {
  return getMatch(match_idx)->getPeptide()->getModifiedAASequence();
} 

int MPSM_Match::getMaxCharge() {

  int max_charge = getCharge(0);

  for (int idx=1;idx < numMatches();idx++) {
    max_charge = max(max_charge, getCharge(idx));
  }

  return max_charge;

}

int MPSM_Match::getCharge(int match_idx) {
  getMatch(match_idx)->getCharge();
}


FLOAT_T MPSM_Match::getScore(SCORER_TYPE_T match_mode) {


  std::map<SCORER_TYPE_T, FLOAT_T>::iterator it = match_scores_.find(match_mode);

  if (it == match_scores_.end()) {
    FLOAT_T score = MPSM_Scorer::score(*this, match_mode);
    setScore(match_mode, score);
    return score;
  } else {
    return it->second;
  }
}

FLOAT_T MPSM_Match::getScoreConst(SCORER_TYPE_T match_mode) const {

  std::map<SCORER_TYPE_T, FLOAT_T>::const_iterator iter = match_scores_.find(match_mode);

  if (iter == match_scores_.end()) {
    carp(CARP_FATAL," Dont have score for %d",(int)match_mode);
  } else {
    return iter->second;
  }
}

void MPSM_Match::setScore(SCORER_TYPE_T match_mode, FLOAT_T score) {

  match_scores_[match_mode] = score;
}


void MPSM_Match::setRank(SCORER_TYPE_T match_mode, int rank) {
  
  match_rank_[match_mode] = rank; 
}

int MPSM_Match::getRank(SCORER_TYPE_T match_mode) {


  std::map<SCORER_TYPE_T, int>::iterator iter = match_rank_.find(match_mode);
  if (iter == match_rank_.end()) {

    carp(CARP_FATAL,"Don't have match rank for %d!",(int)match_mode);

  }
  return iter->second;
}

void MPSM_Match::setBYIonPossible(int b_y_ion_possible) {
  b_y_ion_possible_ = b_y_ion_possible;
}

int MPSM_Match::getBYIonPossible() {

  return b_y_ion_possible_;
}

void MPSM_Match::setBYIonMatched(int b_y_ion_matched) {
  b_y_ion_matched_ = b_y_ion_matched;
}

int MPSM_Match::getBYIonMatched() {
  return b_y_ion_matched_;
}

void MPSM_Match::setBYIonFractionMatched(FLOAT_T b_y_ion_fraction_matched) {
  b_y_ion_fraction_matched_ = b_y_ion_fraction_matched;
}

FLOAT_T MPSM_Match::getBYIonFractionMatched() {
  return b_y_ion_fraction_matched_;
}


void MPSM_Match::getSpectrumNeutralMasses(vector<FLOAT_T>& neutral_masses) {

  neutral_masses.clear();
  Match* first_match = getMatch(0);
  Spectrum* spectrum = first_match->getSpectrum();

  ZStateIndex& zstate_index = getZStateIndex();

  for (int i=0;i<zstate_index.size();i++) {
    FLOAT_T mass = zstate_index[i].getNeutralMass();
    neutral_masses.push_back(mass);
  }
}

void MPSM_Match::getPeptideMasses(vector<FLOAT_T>& peptide_masses) {
  peptide_masses.clear();
  Match* match = NULL;

  for (int idx=0;idx<numMatches();idx++) {
    Match* match = getMatch(idx);
    Peptide* peptide = match->getPeptide();
    FLOAT_T peptide_mass = peptide->getPeptideMass();
    peptide_masses.push_back(peptide_mass);
  }
}

void MPSM_Match::getPeptideModifiedSequences(vector<string>& modified_sequences) {
  modified_sequences.clear();

  Match* match = NULL;

  for (int idx = 0; idx < numMatches(); idx++) {
    Match* match = getMatch(idx);
    Peptide* peptide = match->getPeptide();
    char* seq = peptide->getUnshuffledModifiedSequence();
    string string_seq(seq);
    modified_sequences.push_back(string_seq);
    free(seq);
  }
}

FLOAT_T MPSM_Match::getRTimeMaxDiff() {
  return RetentionPredictor::getStaticRetentionPredictor()->calcMaxDiff(*this);
}

string MPSM_Match::getString() {
  string ans;

  Match* match = NULL;
  for (int idx = 0; idx < numMatches(); idx++) {
    Match* match = getMatch(idx);
    Peptide* peptide = match->getPeptide();
    char* seq = peptide->getModifiedSequenceWithSymbols();
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


void MPSM_Match::setDeltaCN(double delta_cn) {
  //delta_cn_ = delta_cn;
  (void)delta_cn;
}

double MPSM_Match::getDeltaCN() {
  return 0;//delta_cn_;

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

  oss << getMatch(0)->getRank(XCORR);

  for (int idx=1;idx<numMatches();idx++) {
    oss << "," << getMatch(idx)->getRank(XCORR);
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

string MPSM_Match::getProteinIdsString() {

  vector<string> protein_ids;

  for (int idx=0;idx < matches_.size();idx++) {
    Peptide* peptide = matches_[idx]->getPeptide();
    protein_ids.push_back(peptide->getProteinIdsLocations());
    

  }
  return DelimitedFile::splice(protein_ids,';');
}


FLOAT_T MPSM_Match::getXCorrSumDiff() {
  
  FLOAT_T ans = getScore(XCORR);

  for (int idx=0;idx < matches_.size();idx++) {

    ans -= matches_[idx]->getScore(XCORR);
  }

  return ans;

}

FLOAT_T MPSM_Match::getXCorrMaxDiff() {

  FLOAT_T max_score = matches_[0]->getScore(XCORR);

  for (int idx=1;idx < matches_.size();idx++) {
    FLOAT_T current_score = matches_[idx]->getScore(XCORR);
    if (current_score > max_score) {
      max_score = current_score;
    }
  }

  return getScore(XCORR) - max_score;
}

FLOAT_T MPSM_Match::getAreaRatio() {

  if (matches_.size() <= 1) { 
    return 1.0;
  }

  FLOAT_T ans = matches_.at(1)->getZState().getArea() / 
                matches_.at(0)->getZState().getArea();

  return ans;

}

FLOAT_T MPSM_Match::getMZ1RTimeDiff() {

  if (matches_.size() <= 1) {
    return 0.0;
  }
  //cerr << get_match_zstate(matches_.at(1)).getRTime() << endl;
  //cerr << get_match_zstate(matches_.at(0)).getRTime() << endl;
  FLOAT_T ans = matches_.at(1)->getZState().getRTime() -
                matches_.at(0)->getZState().getRTime();

  return ans;
}

void MPSM_Match::getMatchedIons(
  vector<map<Peak*, vector<Ion*> > >& matched_peak_to_ion_vec) {

  FLOAT_T mz_tolerance = 1.0;//get_double_parameter("mz-bin-width");


  int max_charge = min(getMaxCharge(), 
    get_max_ion_charge_parameter("max-ion-charge"));

  vector<IonSeries*> ion_series_vec;
  vector<IonConstraint*> ion_constraints;

  //cerr <<"Creating ion serieses"<<endl;

  for (int seq_idx=0;seq_idx < numMatches();seq_idx++) {

    IonConstraint* ion_constraint = new IonConstraint(
      get_mass_type_parameter("fragment-mass"),
      min(getCharge(seq_idx), max_charge),
      BY_ION,
      false);

    ion_constraints.push_back(ion_constraint);

    char* sequence = getSequence(seq_idx);
    
    IonSeries* ion_series = new IonSeries(
      sequence,
      getCharge(seq_idx),
      ion_constraints.at(seq_idx));

    ion_series->predictIons();

    ion_series_vec.push_back(ion_series);

    //free(sequence);

  }
  map<Peak*, vector<Ion*> >::iterator find_iter;

  vector<int> num_matched_ions;

  //cerr <<"Assigning ions"<<endl;

  for (int seq_idx = 0; seq_idx < numMatches();seq_idx++) {

    //cerr <<"Sequence:"<<peptide_sequences[seq_idx]<<endl;
    
    num_matched_ions.push_back(0);
    //For each ion, find the max_intensity peak with in the the 
    map<Peak*, vector<Ion*> > temp;
    matched_peak_to_ion_vec.push_back(temp);
    map<Peak*, vector<Ion*> > &matched_peak_to_ion = matched_peak_to_ion_vec.at(seq_idx);

    IonSeries* current_ion_series = ion_series_vec.at(seq_idx);

    for (IonIterator ion_iter = current_ion_series -> begin();
      ion_iter != current_ion_series->end();
      ++ion_iter) {

      Ion* ion = *ion_iter;

      Peak* peak = getSpectrum()->getMaxIntensityPeak(ion->getMassZ(), mz_tolerance);
      if (peak != NULL) {  
        num_matched_ions[seq_idx]++;
        find_iter = matched_peak_to_ion.find(peak);

        if (find_iter == matched_peak_to_ion.end()) {
          vector<Ion*> temp_ions;
          matched_peak_to_ion[peak] = temp_ions;
        }

        matched_peak_to_ion[peak].push_back(ion);
      }
    }
    //cerr <<"Peaks Matched:"<<matched_peak_to_ion.size()<<endl;

    //free the ion_series.
    //cerr <<"Free ion series"<<endl;
    delete current_ion_series;
    //and the constraint
    //cerr <<"Free ion constraint"<<endl;
    delete ion_constraints.at(seq_idx);
    

  }

}

void MPSM_Match::getMatchedIonIntensities(
     double& total_ion_intensity,
     double& ion_intensity1,
     double& ion_intensity2,
     double& ion_intensity12) {

  total_ion_intensity = getSpectrum()->getTotalEnergy();


  vector<map<Peak*, vector<Ion*> > > matched_peak_to_ion_vec;
  getMatchedIons(matched_peak_to_ion_vec);

  ion_intensity1 = 0.0;

  map<Peak*, vector<Ion*> >::iterator peak_iter;

  set<Peak*> all_assigned_peaks;
  
  //cerr <<"Summing up intensities"<<endl;
  for (peak_iter = matched_peak_to_ion_vec.at(0).begin();
    peak_iter != matched_peak_to_ion_vec.at(0).end();
    ++peak_iter) {
    
    Peak* peak = peak_iter->first;
    all_assigned_peaks.insert(peak_iter->first);
    ion_intensity1 += peak->getIntensity();

  }
  //cerr <<"ion1:"<<ion_intensity1<<endl;
  ion_intensity2 = 0.0;
  if (matched_peak_to_ion_vec.size() > 1) {

    for (peak_iter = matched_peak_to_ion_vec.at(1).begin();
      peak_iter != matched_peak_to_ion_vec.at(1).end();
      ++peak_iter) {
      Peak* peak = peak_iter->first;
      all_assigned_peaks.insert(peak_iter->first);
      ion_intensity2 += peak->getIntensity();
    }
  }
  //cerr <<"Ion2:"<<ion_intensity2<<endl;
  ion_intensity12 = 0;
  
  for (set<Peak*>::iterator iter = all_assigned_peaks.begin();
    iter != all_assigned_peaks.end();
    ++iter) {

    Peak* peak = *iter;
    ion_intensity12 += peak->getIntensity();
  }
  //cerr <<"Ion12:"<<ion_intensity12<<endl;
}


map<char,char> same;
bool list_generated = false;
void generate_list() {

  same['I'] = 'L';
  same['L'] = 'I';
  same['E'] = 'Q';
  same['Q'] = 'E';
  same['Q'] = 'K';
  same['K'] = 'Q';
  same['D'] = 'N';
  same['N'] = 'D';

  list_generated=true;
}


bool same_char(char &c1, char& c2) {
  if (c1==c2) return true;
  if (!list_generated) {
    generate_list();
  }

  if (same.find(c1) != same.end()) {
    return same[c1] == c2;
  }
  return false;

}


int hamming_dist(Peptide* pep1, Peptide* pep2) {

  //make sure that pep1 is the longest peptide.
  if (pep1->getLength() < pep2->getLength()) {
    return hamming_dist(pep2, pep1);
  }

  char* seq1 = pep1->getSequence();
  char* seq2 = pep2->getSequence();

  //cout<<seq1<<" "<<seq2;
  int dist = pep1->getLength();
  //cout << " " << dist << endl;


  for (int idx=0;idx<pep2->getLength();idx++) {
    if (same_char(seq1[idx],seq2[idx])) {
      //cout<<seq1[idx]<<" is same as "<<seq2[idx]<<endl;
      dist = dist - 1;
    } else {
      //cout<<seq1[idx]<<" is diff as "<<seq2[idx]<<endl;
    }
  }
  //cout<<" "<<dist<<endl;
  free(seq1);
  free(seq2);
  return dist;

}

int hamming_dist(Match* match1, Match* match2) {

  return hamming_dist(match1->getPeptide(), match2->getPeptide());
}


int MPSM_Match::getHammingDist() {

  if (matches_.size() == 1) 
  {
    return matches_[0]->getPeptide()->getLength();
  } else {
    int min_dist = 10000;
    for (int idx1=0;idx1 < matches_.size()-1;idx1++) {

      for (int idx2=idx1+1;idx2 < matches_.size();idx2++) {
        int current = hamming_dist(matches_[idx1], matches_[idx2]);
        if (current < min_dist) {
          min_dist = current;
        }
  
      }
    }


    return min_dist;
  }
}


FLOAT_T MPSM_Match::getTICRatio() {
  
  if (matches_.size() <=1) {
    return 1.0;
  }

  double total_tic;
  double tic_1=1;
  double tic_2=1;
  double tic_12;


  getMatchedIonIntensities(total_tic, tic_1, tic_2, tic_12);

  double ratio = tic_2 / tic_1;
  //cerr <<"TIC ratio:"<<ratio<<endl;
  return ratio;
}


ostream& operator <<(ostream& os, MPSM_Match& match_obj) {

  //cout <<"operator start"<<endl;
  int precision = get_int_parameter("precision");
  os << setprecision(precision);

  //cout <<"getting first_match"<<endl;
  Match* first_match = match_obj.getMatch(0);
  //cout <<"Getting parent"<<endl;
  MPSM_MatchCollection* parent = match_obj.getParent();
  if (parent == NULL) {
    carp(CARP_WARNING, "Match parent is null!FIX");
  }

  //cout <<"Getting spectrum"<<endl;
  Spectrum* spectrum = first_match->getSpectrum();
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

const vector<Match*>& MPSM_Match::getMatches() const {

  return matches_;
}


bool CompareMPSM_MatchVisited::operator() (
  const vector<Match*>& m1, 
  const vector<Match*>& m2
  ) const {

    return compareMPSM_MatchVisited(m1, m2);
}
