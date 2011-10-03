#ifndef MPSM_MATCH_H
#define MPSM_MATCH_H
#include "objects.h"
#include "Match.h"

#include "ZStateIndex.h"

#include <bitset>
#include <map>
#include <vector>


bool MPSM_MatchCompare(Match* m1, Match* m2);

class MPSM_Match {
 protected:

  MPSM_MatchCollection* parent_;

  std::vector<Match*> matches_; //the matches that the multi-match is from.

    
  std::map<SCORER_TYPE_T, FLOAT_T> match_scores_;
  std::map<SCORER_TYPE_T, int> match_rank_;

  int b_y_ion_matched_;
  int b_y_ion_possible_;
  FLOAT_T b_y_ion_fraction_matched_;
  //double delta_cn_;
  void invalidate();

 public:

  MPSM_Match();
  MPSM_Match(Match* match);
  //MPSM_Match(const MPSM_Match& mpsm_match);

  void setParent(MPSM_MatchCollection* parent);
  MPSM_MatchCollection* getParent();

  std::string getSRankString();

  void setDeltaCN(double delta_cn);
  double getDeltaCN();
  double getZScore();
  //void setXCorrRank(int xcorr_rank);
  //int getXCorrRank();

  void setRank(SCORER_TYPE_T match_mode, int rank);
  int getRank(SCORER_TYPE_T match_mode);


  virtual ~MPSM_Match();


  BOOLEAN_T addMatch(Match* match);

  ZStateIndex& getZStateIndex();
  ZStateIndex calcZStateIndex();

  int getHammingDist();


  bool isChargeHomogeneous();

  const vector<Match*>& getMatches() const;

  Match* getMatch(int match_idx) const;
  Match* operator [] (int match_idx);

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

  void setBYIonPossible(int b_y_ion_possible);
  int getBYIonPossible();

  void setBYIonMatched(int b_y_ion_matched);
  int getBYIonMatched();

  void setBYIonFractionMatched(FLOAT_T b_y_ion_fraction_matched);
  FLOAT_T getBYIonFractionMatched();

  void getMatchedIons(
    std::vector<std::map<Peak*, std::vector<Ion*> > >& matched_peak_to_ion_vec);

  void getMatchedIonIntensities(
    double& total_ion_intensity,
    double& ion_intensity1,
    double& ion_intensity2,
    double& ion_intensity12);

  FLOAT_T getTICRatio();

  FLOAT_T getRTimeMaxDiff();


  //Accessor functions for writing the match.
  int getFirstScan();
  std::string getChargeString();
  FLOAT_T getSpectrumPrecursorMZ();
  std::string getPeptideMassString();
  std::string getNeutralMassString();
  std::string getSequenceString();
  int getMatchesPerSpectrum();
  std::string getProteinIdsString();
  FLOAT_T getMZ1RTimeDiff();
  
  


  FLOAT_T getXCorrSumDiff();
  FLOAT_T getXCorrMaxDiff();

  FLOAT_T getAreaRatio();
  bool isDecoy();

    
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
  friend bool compareMPSM_MatchVisited(const std::vector<Match*>& c1, const std::vector<Match*>& c2);



};

struct CompareMPSM_MatchVisited {
  bool operator() (const std::vector<Match*>& m1, const std::vector<Match*>& m2) const;
};

#endif
