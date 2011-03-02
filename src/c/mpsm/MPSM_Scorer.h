#ifndef MPSM_SCORER_H
#define MPSM_SCORER_H
#include "objects.h"

#include "MPSM_Match.h"


class MPSM_Scorer {
  protected:
  /*
    Spectrum* spectrum_;
    int max_charge_;
    SCORER_T* scorer_;
    IonSeries* ion_series_;
    MPSM_Match* current_mpsm_;
    IonConstraint* ion_constraint_;
  */

    
  public:

    static void createIonSeries(MPSM_Match& mpsm_match, IonSeries* ion_series);
    static void createIonSeries(std::vector<string>& sequences, std::vector<int>& charges, IonSeries* ion_series);

    MPSM_Scorer();

    virtual ~MPSM_Scorer();

    static FLOAT_T score(
      MPSM_Match& mpsm_match, 
      SCORER_TYPE_T match_mode
      );

};

#endif
