#ifndef XLINKMATCH_H_
#define XLINKMATCH_H_


#include "crux-utils.h"

#include "XLink.h"
#include "XLinkMatchCollection.h"

#include <string>


enum _xlinkmatch_type  
  {LINEAR_CANDIDATE, 
   SELFLOOP_CANDIDATE, 
   XLINK_CANDIDATE};

typedef enum _xlinkmatch_type XLINKMATCH_TYPE_T; 

class XLinkMatch : public Match {

 protected:
  XLinkMatchCollection* parent_;
  FLOAT_T xcorr_;
  int xcorr_rank_;
  FLOAT_T sp_;
  int sp_rank_;

  int by_ions_matched_;
  int by_ions_total_;

  FLOAT_T pvalue_;
  bool mass_calculated_[NUMBER_MASS_TYPES];
  FLOAT_T mass_[NUMBER_MASS_TYPES];
 public:
  XLinkMatch();

  virtual ~XLinkMatch();

  virtual XLINKMATCH_TYPE_T getCandidateType() = 0;

  virtual std::string getSequenceString() = 0;
  virtual FLOAT_T calcMass(MASS_TYPE_T mass_type) = 0;
  
  FLOAT_T getMass(MASS_TYPE_T mass_type);

  virtual XLinkMatch* shuffle() = 0;

  virtual void predictIons(IonSeries* ion_series, int charge)=0;
  virtual std::string getIonSequence(Ion* ion)=0;
  virtual PEPTIDE_T* getPeptide(int peptide_idx)=0;


  void computeWeibullPvalue(
    FLOAT_T shift,
    FLOAT_T eta,
    FLOAT_T beta);

  
  void setXCorr(FLOAT_T xcorr);
  FLOAT_T getXCorr();

  void setXCorrRank(int xcorr_rank);
  int getXCorrRank();

  void setSP(FLOAT_T sp);
  FLOAT_T getSP();

  void setSPRank(int sp_rank);
  int getSPRank();

  void setBYIonsMatched(int by_ions_matched);
  int getBYIonsMatched();

  void setBYIonsTotal(int by_ions_total);
  int getBYIonsTotal();

  virtual int getNumMissedCleavages() = 0;

  virtual bool isModified() = 0;


  static std::string getResultHeader();
  std::string getResultString();
  FLOAT_T getPPMError();
  std::string getProteinIdString(int peptide_idx);
  void setParent(XLinkMatchCollection* parent);

  


};

#endif
