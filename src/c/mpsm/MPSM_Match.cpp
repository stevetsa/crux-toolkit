#include "MPSM_Match.h"
#include "MPSM_Scorer.h"

using namespace std;

MPSM_Match::MPSM_Match() {
  invalidate();
}

MPSM_Match::MPSM_Match(MATCH_T* match) {
  addMatch(match);
}

MPSM_Match::MPSM_Match(const MPSM_Match& mpsm_match) {
  //TODO: copy match over.
}


MPSM_Match::~MPSM_Match() {
  matches_.clear();
}

BOOLEAN_T MPSM_Match::addMatch(MATCH_T* match) {

  //check to make sure match doesn't already exist here.
  for (int idx=0;idx<matches_.size();idx++) {
    if (matches_[idx] == match)
      return FALSE;
  }
  matches_.push_back(match);
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

void MPSM_Match::invalidate() {
  for (int idx=0;idx < _SCORE_TYPE_NUM;idx++) {
    match_scores_valid_[idx] = FALSE;
  }
  charge_valid_ = FALSE;
}


ChargeIndex& MPSM_Match::getChargeIndex() {
  if (!charge_valid_) {
    //TODO: validate the charge_index.

  }

  return charge_index_;
}

MATCH_T* MPSM_Match::getMatch(int match_idx) {
  return matches_[match_idx];
}

int MPSM_Match::numMatches() {
  return matches_.size();
}


FLOAT_T MPSM_Match::getScore(SCORER_TYPE_T match_mode) {

  if (!match_scores_valid_[match_mode]) {
    FLOAT_T score = MPSM_Scorer::scoreMPSM(*this, match_mode);
    setScore(match_mode, score);
    return score;
  }
  return match_scores_[match_mode];
}

void MPSM_Match::setScore(SCORER_TYPE_T match_mode, FLOAT_T score) {

  match_scores_[match_mode] = score;
  match_scores_valid_[match_mode] = TRUE;
}
