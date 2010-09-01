#ifndef MATCHCANDIDATE_H_
#define MATCHCANDIDATE_H_


#include "crux-utils.h"

#include "XLink.h"
#include "MatchCandidateVector.h"

#include <string>


enum _matchcandidate_type  
  {LINEAR_CANDIDATE, 
   SELFLOOP_CANDIDATE, 
   XLINK_CANDIDATE};

typedef enum _matchcandidate_type MATCHCANDIDATE_TYPE_T; 

class MatchCandidate {

 protected:
  MatchCandidateVector* parent_;
  FLOAT_T xcorr_;
  int xcorr_rank_;
  FLOAT_T pvalue_;
 public:

  virtual ~MatchCandidate();

  virtual MATCHCANDIDATE_TYPE_T getCandidateType() = 0;

  virtual std::string getSequenceString() = 0;
  virtual FLOAT_T getMass() = 0;

  virtual MatchCandidate* shuffle() = 0;

  virtual void predictIons(ION_SERIES_T* ion_series, int charge)=0;
  virtual std::string getIonSequence(ION_T* ion)=0;

  void computeWeibullPvalue(
    FLOAT_T shift,
    FLOAT_T eta,
    FLOAT_T beta);

  
  void setXCorr(FLOAT_T xcorr);
  FLOAT_T getXCorr();

  static std::string getResultHeader();
  std::string getResultString();

  void setParent(MatchCandidateVector* parent);


};

#endif
