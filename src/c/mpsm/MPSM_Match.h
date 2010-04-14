#ifndef MPSM_MATCH_H
#define MPSM_MATCH_H
#include "objects.h"
#include "match.h"

#include "ChargeIndex.h"

#include <vector>


bool MPSM_MatchCompare(MATCH_T* m1, MATCH_T* m2);

class MPSM_Match {
  protected:

    MPSM_MatchCollection* parent_;

    std::vector<MATCH_T*> matches_; //the matches that the multi-match is from.
    FLOAT_T match_scores_[_SCORE_TYPE_NUM];
    BOOLEAN_T match_scores_valid_[_SCORE_TYPE_NUM];
    
    BOOLEAN_T charge_valid_;

    ChargeIndex charge_index_;

    FLOAT_T rtime_max_diff_;

    FLOAT_T delta_cn_;

    std::vector<BOOLEAN_T> has_rtime;
    std::vector<FLOAT_T> rtimes;
    void invalidate();

  public:

    MPSM_Match();
    MPSM_Match(MATCH_T* match);
    //MPSM_Match(const MPSM_Match& mpsm_match);

    void setParent(MPSM_MatchCollection* parent);
    MPSM_MatchCollection* getParent();



    virtual ~MPSM_Match();


    BOOLEAN_T addMatch(MATCH_T* match);

    BOOLEAN_T isDecoy();

    BOOLEAN_T hasRTime(int match_idx);
    FLOAT_T getRTime(int match_idx);
    void setRTime(int match_idx, FLOAT_T rtime);

    void setRTimeMaxDiff(FLOAT_T rtime_max_diff);
    FLOAT_T getRTimeMaxDiff();

    void setDeltaCN(FLOAT_T delta_cn);


    ChargeIndex& getChargeIndex();

    MATCH_T* getMatch(int match_idx);
    MATCH_T* operator [] (int match_idx);


    int numMatches();

    
    FLOAT_T getScore(SCORER_TYPE_T match_mode);
    void setScore(SCORER_TYPE_T match_mode, FLOAT_T score);

    void getSpectrumNeutralMasses(std::vector<FLOAT_T>& neutral_masses);
    void getPeptideMasses(std::vector<FLOAT_T>& peptide_masses);
    void getPeptideModifiedSequences(std::vector<std::string>& modified_sequences);

    static void sortMatches(std::vector<MPSM_Match>& matches, SCORER_TYPE_T match_mode);
    
    bool operator == (MPSM_Match& match_obj);
    bool operator < (const MPSM_Match& match_obj) const;

    friend std::ostream& operator <<(std::ostream& os, MPSM_Match& match_obj);

    friend bool compareMPSM_Match(const MPSM_Match& c1, const MPSM_Match& c2);


};

#endif
