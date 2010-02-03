#ifndef MPSM_MATCH_H
#define MPSM_MATCH_H
#include "objects.h"
#include "match.h"

#include "ChargeIndex.h"

#include <vector>


class MPSM_Match {
  protected:
    std::vector<MATCH_T*> matches_; //the matches that the multi-match is from.
    FLOAT_T match_scores_[_SCORE_TYPE_NUM];
    BOOLEAN_T match_scores_valid_[_SCORE_TYPE_NUM];

    BOOLEAN_T charge_valid_;

    ChargeIndex charge_index_;

    void invalidate();

  public:

    MPSM_Match();
    MPSM_Match(MATCH_T* match);
    MPSM_Match(const MPSM_Match& mpsm_match);

    virtual ~MPSM_Match();

    BOOLEAN_T addMatch(MATCH_T* match);

    BOOLEAN_T isDecoy();

    ChargeIndex& getChargeIndex();

    MATCH_T* getMatch(int match_idx);
    int numMatches();

    
    FLOAT_T getScore(SCORER_TYPE_T match_mode);
    void setScore(SCORER_TYPE_T match_mode, FLOAT_T score);

};

#endif
