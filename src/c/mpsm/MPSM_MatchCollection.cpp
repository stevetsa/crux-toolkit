#include "MPSM_MatchCollection.h"
#include "RetentionPredictor.h"
#include "MPSM_Scorer.h"

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
  if (sorted_ && sort_mode_ == match_mode) return;
  //sort by score.
  //cerr <<"Scoring"<<endl;

  if (matches_.size() == 0) {
    sort_mode_ = match_mode;
    sorted_ = true;
    return;
  }
  MPSM_Scorer scorer(matches_[0].getSpectrum(), matches_[0].getMaxCharge(), match_mode);
  for (int i=0;i<numMatches();i++) {
    FLOAT_T score = scorer.calcScore(matches_[i]);
    matches_[i].setScore(match_mode, score);
  }
  //cerr <<"Sorting"<<endl;
  MPSM_Match::sortMatches(matches_, match_mode);
  sorted_ = true;
  sort_mode_ = match_mode;
  //cerr <<"Done"<<endl;
}

void MPSM_MatchCollection::calcDeltaCN() {

  if (numMatches() == 0) {
    xcorr_2_ = 0;
    return;
  }

  if (numMatches() == 1) {
    xcorr_2_ = 0;
  } else {

    sortByScore(XCORR);  
    xcorr_2_ = getMatch(1).getScore(XCORR);
  }
}

void MPSM_MatchCollection::calcRanks(SCORER_TYPE_T scorer_type) {
  if (numMatches() == 0) {
    return;
  }
/*
  else if (getMatch(0).numMatches() == 1) {
    //all single peptide match ranks are already calculated.
    for (int idx=0;idx < numMatches();idx++) {
      int rank = get_match_rank(getMatch(idx).getMatch(0), scorer_type);
      getMatch(idx).setRank(scorer_type, rank);
    }
  }*/ else {
    int top_match = 0;

    if (scorer_type == XCORR) {
      top_match = get_int_parameter("mpsm-top-match");
    }  
    sortByScore(scorer_type);
    
    int cur_rank = 1;
    FLOAT_T cur_score = getMatch(0).getScore(scorer_type);
    getMatch(0).setRank(scorer_type, cur_rank);

    for (int match_idx = 1; match_idx < numMatches();match_idx++) {
      FLOAT_T this_score = getMatch(match_idx).getScore(scorer_type);
      if (this_score != cur_score) {
        cur_score = this_score;
        cur_rank++;
      }
      getMatch(match_idx).setRank(scorer_type, cur_rank);
      if ((scorer_type == XCORR) && (cur_rank > top_match)) break;
    }

  }

}



double MPSM_MatchCollection::calcDeltaCNMatch(double xcorr) {
  if (xcorr > 0 && numMatches() > 1) {
    return (xcorr - xcorr_2_) / xcorr;
  } else {
    return 0;
  }
}


void MPSM_MatchCollection::calcZParameters(double& mean, double& std) {
  mean = 0;
  std = 1;

  if (numMatches() < 2) {
    return;
  }

  for (int idx=0;idx < numMatches();idx++) {
    mean += getMatch(idx).getScore(XCORR);
  }

  mean = mean / (double)numMatches();

  std = 0;
  for (int idx=0;idx < numMatches();idx++) {
    double temp = getMatch(idx).getScore(XCORR) - mean;
    std += temp * temp;
  }

  std = sqrt(1.0 / (double)(numMatches() - 1) * std);

  setZParameters(mean, std);

}

void MPSM_MatchCollection::setZParameters(double mean, double std) {
  zscore_mean_ = mean;
  zscore_std_ = std;
}

double MPSM_MatchCollection::calcZScore(double xcorr) {

    return (xcorr - zscore_mean_) / zscore_std_;
}

double MPSM_MatchCollection::getSpectrumRTime() {
  return 0.0;
  //return get_spectrum_rtime(get_match_spectrum(getMatch(0)[0]));
}

RetentionPredictor* rtime_predictor_ = NULL;

double MPSM_MatchCollection::getPredictedRTime(MPSM_Match& match) {
  double ans = 0.0;

  if (match.numMatches() == 1) {
    if (rtime_predictor_ == NULL) {
      rtime_predictor_ = RetentionPredictor::createRetentionPredictor();
    }
    //cerr<<"Predicting retention time"<<endl;
    ans = rtime_predictor_->predictRTime(match[0]);
  }
  return ans;

  
}

bool MPSM_MatchCollection::visited(
  MPSM_Match& match
  ) {

  pair<set<MPSM_Match>::iterator,bool> ret;
  ret = visited_.insert(match);
  
  return !ret.second;
}



ostream& operator<<(ostream& os, MPSM_MatchCollection& collection_obj) {

  for (int idx=0;idx < collection_obj.numMatches();idx++) {
    os << collection_obj.getMatch(idx) << endl;
  }

  return os;
}
