#ifndef MPSM_MATCH_H
#define MPSM_MATCH_H
#include "objects.h"
#include "match.h"

#include "ZStateIndex.h"

#include <bitset>
#include <map>
#include <vector>


bool MPSM_MatchCompare(MATCH_T* m1, MATCH_T* m2);

class MPSM_Match {
 protected:

  MPSM_MatchCollection* parent_;

  std::vector<MATCH_T*> matches_; //the matches that the multi-match is from.

    
  FLOAT_T match_scores_[NUMBER_SCORER_TYPES];
  std::map<SCORER_TYPE_T, int> match_rank_;
  std::bitset<NUMBER_SCORER_TYPES> have_match_score_;
  std::bitset<NUMBER_SCORER_TYPES> have_match_rank_;

  bool zstate_valid_;

  ZStateIndex zstate_index_;

    
  FLOAT_T rtime_max_diff_;

  /*
  std::vector<bool> has_rtime;
  std::vector<FLOAT_T> rtimes;
  */
  void invalidate();

 public:

  MPSM_Match();
  MPSM_Match(MATCH_T* match);
  //MPSM_Match(const MPSM_Match& mpsm_match);

  void setParent(MPSM_MatchCollection* parent);
  MPSM_MatchCollection* getParent();

  std::string getSRankString();

  double getDeltaCN();
  double getZScore();
  //void setXCorrRank(int xcorr_rank);
  //int getXCorrRank();

  void setRank(SCORER_TYPE_T match_mode, int rank);
  int getRank(SCORER_TYPE_T match_mode);


  virtual ~MPSM_Match();


  BOOLEAN_T addMatch(MATCH_T* match);

  BOOLEAN_T isDecoy();

  BOOLEAN_T hasRTime(int match_idx);
  FLOAT_T getRTime(int match_idx);
  void setRTime(int match_idx, FLOAT_T rtime);

  void setRTimeMaxDiff(FLOAT_T rtime_max_diff);
  FLOAT_T getRTimeMaxDiff();

  FLOAT_T getRTimeAveDiff();


  ZStateIndex& getZStateIndex();

  bool isChargeHomogeneous();


  MATCH_T* getMatch(int match_idx) const;
  MATCH_T* operator [] (int match_idx);

  Spectrum* getSpectrum();
  char* getSequence(int match_idx);
  int getCharge(int match_idx);

  MODIFIED_AA_T* getModifiedSequence(int match_idx);

  int getMaxCharge();

  double getSpectrumRTime();
  double getPredictedRTime();
  int numMatches() const;


  std::string getString();

  void predictIons(IonSeries* ion_series);

  //Accessor functions for writing the match.
  int getFirstScan();
  std::string getChargeString();
  FLOAT_T getSpectrumPrecursorMZ();
  std::string getPeptideMassString();
  std::string getNeutralMassString();
  std::string getSequenceString();
  int getMatchesPerSpectrum();

  FLOAT_T getXCorrSumDiff();
  FLOAT_T getXCorrMaxDiff();

  FLOAT_T getAreaRatio();


    
  FLOAT_T getScore(SCORER_TYPE_T match_mode);
  FLOAT_T getScoreConst(SCORER_TYPE_T match_mode) const;
  void setScore(SCORER_TYPE_T match_mode, FLOAT_T score);
  void getSpectrumNeutralMasses(std::vector<FLOAT_T>& neutral_masses);
  void getPeptideMasses(std::vector<FLOAT_T>& peptide_masses);
  void getPeptideModifiedSequences(std::vector<std::string>& modified_sequences);
  static void sortMatches(std::vector<MPSM_Match>& matches, SCORER_TYPE_T match_mode);
    
  bool operator == (MPSM_Match& match_obj);
  bool operator < (const MPSM_Match& match_obj) const;
  friend std::ostream& operator <<(std::ostream& os, MPSM_Match& match_obj);
  friend bool compareMPSM_Match(const MPSM_Match& c1, const MPSM_Match& c2);
  friend bool compareMPSM_MatchVisited(const MPSM_Match& c1, const MPSM_Match& c2);



};

struct CompareMPSM_MatchVisited {
  bool operator() (const MPSM_Match& m1, const MPSM_Match& m2) const;
};

#endif
