#ifndef MATCHSEARCH_H
#define MATCHSEARCH_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"
#include "OutputFiles.h"

#include <string>

class MatchSearch: public CruxApplication {

 protected:
  /* Private functions */
  int search_pep_mods(
    MATCH_COLLECTION_T* match_collection, ///< store PSMs here
    bool is_decoy,   ///< generate decoy peptides from index/db
    INDEX_T* index,       ///< index to use for generating peptides
    DATABASE_T* database, ///< db to use for generating peptides
    Spectrum* spectrum, ///< spectrum to search
    int charge,           ///< seach spectrum at this charge state
    PEPTIDE_MOD_T** pep_mod_list, ///< list of peptide mods to apply
    int num_peptide_mods, ///< how many p_mods to use from the list
    bool store_scores///< keep all scores for p-value estimation
  );

  void add_decoy_scores(
    MATCH_COLLECTION_T* target_psms, ///< add scores to these matches
    Spectrum* spectrum, ///<
    int charge, ///< 
    INDEX_T* index, ///< search this index if not null
    DATABASE_T* database, ///< search this database if not null
    PEPTIDE_MOD_T** peptitde_mods, ///< list of peptide mods to search
    int num_peptide_mods ///< number of mods in the above array
  );

  bool is_search_complete(MATCH_COLLECTION_T* matches, 
                             int mods_per_peptide);

  void print_spectrum_matches(
    OutputFiles& output_files,       
    MATCH_COLLECTION_T* target_psms, 
    MATCH_COLLECTION_T** decoy_psms,
    int num_decoy_collections,
    Spectrum* spectrum,             
    bool combine_target_decoy,
    int num_decoy_files
  );

 public:

  MatchSearch();
  ~MatchSearch();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();
  virtual std::string getFileString();

};


#endif
