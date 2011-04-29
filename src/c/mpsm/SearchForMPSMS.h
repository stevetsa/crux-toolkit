/**
 * \file SearchForMPSMS.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Object for running search-for-xlinks
 *****************************************************************************/

#ifndef SEARCHFORMPSMS_H
#define SEARCHFORMPSMS_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"

#include "MPSM_ZStateMap.h"
#include "RetentionPredictor.h"

#include <string>

class SearchForMPSMS: public CruxApplication {

 protected:
  bool rtime_threshold_;
  double rtime_all2_threshold_;
  double rtime_all3_threshold_;
  double rtime_default_threshold_;

  RetentionPredictor* rtime_predictor_;

  void search(
    MPSM_ZStateMap& charge_spsm_map, 
    MPSM_ZStateMap& charge_mpsm_map
  );

  int searchPepMods(
    MATCH_COLLECTION_T* match_collection, ///< store PSMs here
    BOOLEAN_T is_decoy,   ///< generate decoy peptides from index/db
    INDEX_T* index,       ///< index to use for generating peptides
    DATABASE_T* database, ///< db to use for generating peptides
    Spectrum* spectrum,         ///< spectrum to search
    SpectrumZState& zstate,       ///< seach spectrum at this charge state
    PEPTIDE_MOD_T** peptide_mods, ///< list of peptide mods to apply
    int num_peptide_mods, ///< how many p_mods to use from the list
    BOOLEAN_T store_scores///< keep all scores for p-value estimation
  );

  bool isSearchComplete(
    MATCH_COLLECTION_T* matches, 
    int mods_per_peptide
  );

  bool extendChargeMap(
    MPSM_ZStateMap& spsm_map,
    MPSM_ZStateMap& current_mpsm_map,
    int mpsm_level
  );

  bool extendMatch(
    MPSM_Match& orig_mpsm,
    MPSM_MatchCollection& spsm_matches,
    MPSM_ZStateMap& new_mpsm_matches,
    int match_collection_idx);

  
  bool passRTimeThreshold(
    MPSM_Match& match
  );

 public:

  SearchForMPSMS();
  ~SearchForMPSMS();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();

  /**
   * \returns the file stem of the application, default getName.
   */ 
  virtual std::string getFileStem();

  /**
   * \returns the enum of the application, default MISC_COMMAND
   */
  virtual COMMAND_T getCommand();

  /**
   * \returns whether the application needs the output directory or not. (default false).
   */
  virtual bool needsOutputDirectory();



};


#endif
