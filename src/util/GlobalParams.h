#include "objects.h"
#include <string>
#include <vector>


class GlobalParams {

 protected:
  static MASS_TYPE_T isotopic_mass_;
  static int missed_cleavages_;
  static int max_aas_modified_;
  static FLOAT_T min_mass_;
  static FLOAT_T max_mass_;
  static WINDOW_TYPE_T precursor_window_type_;
  static FLOAT_T precursor_window_;
  static int min_length_;
  static int max_length_;
  static std::string xlink_prevents_cleavage_;
  static std::string max_ion_charge_;
  static ION_TYPE_T primary_ions_;
  static MASS_TYPE_T fragment_mass_;
  static bool precursor_ions_;
  static ENZYME_T enzyme_;
  static DIGEST_T digestion_;
  static FLOAT_T remove_precursor_tolerance_;
  static OBSERVED_PREPROCESS_STEP_T stop_after_;
  static bool xlink_include_inter_;
  static bool xlink_include_intra_;
  static bool xlink_include_inter_intra_;
  static int max_xlink_mods_;
  static int mod_precision_;
  static int xlink_top_n_;
  static std::vector<int> isotope_windows_;

 public:
  static void set();
  
  static const MASS_TYPE_T& getIsotopicMass();
  static const int& getMissedCleavages();
  static const int& getMaxAasModified();
  static const FLOAT_T& getMinMass();
  static const FLOAT_T& getMaxMass();
  static const WINDOW_TYPE_T& getPrecursorWindowType();
  static const FLOAT_T& getPrecursorWindow();
  static const int& getMinLength();
  static const int& getMaxLength();
  static const std::string& getXLinkPreventsCleavage();
  static const std::string& getMaxIonCharge();
  static const ION_TYPE_T& getPrimaryIons();
  static const MASS_TYPE_T& getFragmentMass();
  static const bool& getPrecursorIons();
  static const ENZYME_T& getEnzyme();
  static const DIGEST_T& getDigestion();
  static const FLOAT_T& getRemovePrecursorTolerance();
  static const OBSERVED_PREPROCESS_STEP_T& getStopAfter();
  static const bool& getXlinkIncludeInter();
  static const bool& getXLinkIncludeIntra();
  static const bool& getXLinkIncludeInterIntra();
  static const int& getMaxXLinkMods();
  static const int& getModPrecision();
  static const int& getXLinkTopN();
  static const std::vector<int>& getIsotopeWindows();

};
