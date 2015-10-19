#include "GlobalParams.h"
#include "parameter.h"
#include "util/StringUtils.h"

using namespace std;


MASS_TYPE_T GlobalParams::isotopic_mass_;
int GlobalParams::missed_cleavages_;
int GlobalParams::max_aas_modified_;
FLOAT_T GlobalParams::min_mass_;
FLOAT_T GlobalParams::max_mass_;
WINDOW_TYPE_T GlobalParams::precursor_window_type_;
FLOAT_T GlobalParams::precursor_window_;
int GlobalParams::min_length_;
int GlobalParams::max_length_;
string GlobalParams::xlink_prevents_cleavage_;
string GlobalParams::max_ion_charge_;
ION_TYPE_T GlobalParams::primary_ions_;
MASS_TYPE_T GlobalParams::fragment_mass_;
bool GlobalParams::precursor_ions_;
ENZYME_T GlobalParams::enzyme_;
DIGEST_T GlobalParams::digestion_;
FLOAT_T GlobalParams::remove_precursor_tolerance_;
OBSERVED_PREPROCESS_STEP_T GlobalParams::stop_after_;
bool GlobalParams::xlink_include_inter_;
bool GlobalParams::xlink_include_intra_;
bool GlobalParams::xlink_include_inter_intra_;
int GlobalParams::max_xlink_mods_;
int GlobalParams::mod_precision_;
int GlobalParams::xlink_top_n_;
vector<int> GlobalParams::isotope_windows_;


void GlobalParams::set() {
  bool monoisotopic_precursor = 
    get_boolean_parameter("monoisotopic-precursor");
  if (monoisotopic_precursor) {
    isotopic_mass_ = MONO;
  } else {
    isotopic_mass_ = AVERAGE;
  }

  missed_cleavages_ = get_int_parameter("missed-cleavages");
  max_aas_modified_ = get_int_parameter("max-aas-modified");
  min_mass_ = get_double_parameter("min-mass");
  max_mass_ = get_double_parameter("max-mass");
  precursor_window_type_ = string_to_window_type(get_string_parameter("precursor-window-type"));
  precursor_window_ = get_double_parameter("precursor-window");
  min_length_ = get_int_parameter("min-length");
  max_length_ = get_int_parameter("max-length");
  xlink_prevents_cleavage_ = get_string_parameter("xlink-prevents-cleavage");
  max_ion_charge_ = get_string_parameter("max-ion-charge");
  string_to_ion_type(get_string_parameter("primary-ions"), &primary_ions_);
  fragment_mass_ = get_mass_type_parameter("fragment-mass");
  precursor_ions_ = get_boolean_parameter("precursor-ions");
  enzyme_ = get_enzyme_type_parameter("enzyme");
  digestion_ = get_digest_type_parameter("digestion");
  remove_precursor_tolerance_ = get_double_parameter("remove-precursor-tolerance");
  stop_after_ = string_to_observed_preprocess_step(get_string_parameter("stop-after"));
  xlink_include_inter_ = get_boolean_parameter("xlink-include-inter");
  xlink_include_intra_ = get_boolean_parameter("xlink-include-intra");
  xlink_include_inter_intra_ = get_boolean_parameter("xlink-include-inter-intra");
  max_xlink_mods_ = get_int_parameter("max-xlink-mods");
  mod_precision_ = get_int_parameter("mod-precision");
  xlink_top_n_ = get_int_parameter("xlink-top-n");
  isotope_windows_ = StringUtils::Split<int>(get_string_parameter("isotope-windows"), ',');
}

const MASS_TYPE_T& GlobalParams::getIsotopicMass() {
  return isotopic_mass_;
}

const int& GlobalParams::getMissedCleavages() {
  return missed_cleavages_;
}

const int& GlobalParams::getMaxAasModified() {
  return max_aas_modified_;
}

const FLOAT_T& GlobalParams::getMinMass() {
  return min_mass_;
}

const FLOAT_T& GlobalParams::getMaxMass() {
  return max_mass_;
}

const WINDOW_TYPE_T& GlobalParams::getPrecursorWindowType() {
  return precursor_window_type_;
}

const FLOAT_T& GlobalParams::getPrecursorWindow() {
  return precursor_window_;
}

const int& GlobalParams::getMinLength() {
  return min_length_;
}

const int& GlobalParams::getMaxLength() {
  return max_length_;
}

const string& GlobalParams::getXLinkPreventsCleavage() {
  return xlink_prevents_cleavage_;
}

const string& GlobalParams::getMaxIonCharge() {
  return max_ion_charge_;
}

const ION_TYPE_T& GlobalParams::getPrimaryIons() {
  return primary_ions_;
}

const MASS_TYPE_T& GlobalParams::getFragmentMass() {
  return fragment_mass_;
}

const bool& GlobalParams::getPrecursorIons() {
  return precursor_ions_;
}

const ENZYME_T& GlobalParams::getEnzyme() {
  return enzyme_;
}

const DIGEST_T& GlobalParams::getDigestion() {
  return digestion_;
}

const FLOAT_T& GlobalParams::getRemovePrecursorTolerance() {
  return remove_precursor_tolerance_;
}

const OBSERVED_PREPROCESS_STEP_T& GlobalParams::getStopAfter() {
  return stop_after_;
}

const bool& GlobalParams::getXlinkIncludeInter() {
  return xlink_include_inter_;
}
  
const bool& GlobalParams::getXLinkIncludeIntra() {
  return xlink_include_intra_;
}

const bool& GlobalParams::getXLinkIncludeInterIntra() {
  return xlink_include_inter_intra_;
}
  
const int& GlobalParams::getMaxXLinkMods() {
  return max_xlink_mods_;
}

const int& GlobalParams::getModPrecision() {
  return mod_precision_;
}

const int& GlobalParams::getXLinkTopN() {
  return xlink_top_n_;
}

const vector<int>& GlobalParams::getIsotopeWindows() {
  return isotope_windows_;
}
