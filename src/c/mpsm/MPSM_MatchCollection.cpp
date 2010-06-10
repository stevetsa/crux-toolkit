#include "MPSM_MatchCollection.h"

#include <iostream>

using namespace std;


MPSM_MatchCollection::MPSM_MatchCollection() {
  spsm_matches_ = NULL;
  sorted_ = false;

}

MPSM_MatchCollection::MPSM_MatchCollection(MATCH_COLLECTION_T* spsm_matches) {

  spsm_matches_ = spsm_matches;
  MATCH_ITERATOR_T *match_iter =
    new_match_iterator(spsm_matches, XCORR, FALSE);

  while(match_iterator_has_next(match_iter)) {
    MATCH_T* match = match_iterator_next(match_iter);
    MPSM_Match mpsm_match(match);
    addMatch(mpsm_match);
    /*
    if (mpsm_match.getParent() == NULL) {
      carp(CARP_FATAL,"Why is this null?");
    }
    if (getMatch(numMatches()-1).getParent() == NULL) {
      carp(CARP_FATAL,"Why is last null?");
    }
    */
  }
  free_match_iterator(match_iter);
  sorted_ = false;

  //cout <<"Num matches added:"<<matches_.size()<<endl;
  //cout <<"Matches total:"<<get_match_collection_match_total(spsm_matches)<<endl;

}

MPSM_MatchCollection::~MPSM_MatchCollection() {
  matches_.clear();
}

void MPSM_MatchCollection::free() {
  if (spsm_matches_ != NULL) {
    free_match_collection(spsm_matches_);
    spsm_matches_ = NULL;
  }
}

bool MPSM_MatchCollection::addMatch(MPSM_Match& match) {
  match.setParent(this);
  matches_.push_back(match);
  sorted_ = false;
  return TRUE;

}


MPSM_Match& MPSM_MatchCollection::getMatch(int idx) {
  return matches_[idx];
}
  
MPSM_Match& MPSM_MatchCollection::operator [] (int idx) {
  return getMatch(idx);
}
  
int MPSM_MatchCollection::numMatches() {
  return matches_.size();
}

void MPSM_MatchCollection::sortByScore(SCORER_TYPE_T match_mode) {
  if (sorted_) return;
  //sort by score.

  for (int i=0;i<numMatches();i++) {
    matches_[i].getScore(match_mode);
  }

  MPSM_Match::sortMatches(matches_, match_mode);
  sorted_ = true;
}

void MPSM_MatchCollection::calcDeltaCN() {

  if (numMatches() == 0) {
    return;
  }

  if (numMatches() == 1) {

    getMatch(0).setDeltaCN(0.0);
  } else {

    sortByScore(XCORR);
  
    FLOAT_T second_best = getMatch(1).getScore(XCORR);

    for (int idx = 0;idx < numMatches();idx++) {

      FLOAT_T delta_cn = (getMatch(idx).getScore(XCORR) - second_best) / 
        second_best;
      getMatch(idx).setDeltaCN(delta_cn);
    }
  }
}

void MPSM_MatchCollection::calcZScores() {
  if (numMatches() == 0) {
    return;
  }

  if (numMatches() < 2) {
    for (int idx=0;idx < numMatches();idx++) {
      getMatch(idx).setZScore(0);
    }
  } else {
    double sum = 0;

    for (int idx = 0;idx < numMatches();idx++) {
      sum += getMatch(idx).getScore(XCORR);
    }

    double avg = sum / (double)numMatches();

    double std = 0;

    for (int idx=0;idx < numMatches();idx++) {
      double temp = getMatch(idx).getScore(XCORR) - avg;
      std += temp * temp;

    }

    std = sqrt(1.0 / (double)(numMatches() - 1) * std);

    for (int idx=0;idx < numMatches();idx++) {

      double zscore = (getMatch(idx).getScore(XCORR) - avg) / std;
      getMatch(idx).setZScore(zscore);


    }


  }


}


ostream& operator<<(ostream& os, MPSM_MatchCollection& collection_obj) {

  for (int idx=0;idx < collection_obj.numMatches();idx++) {
    os << collection_obj.getMatch(idx) << endl;
  }

  return os;
}
