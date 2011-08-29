#ifndef MPSM_SCORER_H
#define MPSM_SCORER_H
#include "objects.h"

#include "MPSM_Match.h"


class MPSM_Scorer {
  protected:
  
    Spectrum* spectrum_;
    int max_charge_;
    Scorer* scorer_;
    SCORER_TYPE_T scorer_type_;
    IonConstraint* ion_constraint_;
  /*
    IonSeries* ion_series_;
    MPSM_Match* current_mpsm_;
    IonConstraint* ion_constraint_;
  */
    //static vector<vector< IonConstraint*> > ion_constraints_;
    static IonConstraint* getIonConstraint(SCORER_TYPE_T scorer_type, int charge);
    static IonSeries* getIonSeries(IonConstraint* ion_constraint, int charge);
  public:

    static void createIonSeries(MPSM_Match& mpsm_match, IonSeries* ion_series);
    static void createIonSeries(std::vector<string>& sequences, std::vector<int>& charges, IonSeries* ion_series);

    MPSM_Scorer();

    MPSM_Scorer(
        Spectrum* spectrum,
        int max_charge,
        SCORER_TYPE_T scorer_type);
        

    FLOAT_T calcScore(MPSM_Match& mpsm_match);

    virtual ~MPSM_Scorer();

    static FLOAT_T score(
      MPSM_Match& mpsm_match, 
      SCORER_TYPE_T match_mode
      );

};

#endif
