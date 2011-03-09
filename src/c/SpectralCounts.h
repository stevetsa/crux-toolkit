/**
 * \file SpectralCounts.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Object for printing the crux version number.
 *****************************************************************************/
#ifndef SPECRAL_COUNTS_H
#define SPECRAL_COUNTS_H

//#include "CruxApplication.h" // restore this line when merged into trunk
#include "spectral-counts.h"
#include <string>
#include <vector>
#include <map>
#include <set>
#include "utils.h"
#include "objects.h"



//class SpectralCounts: public CruxApplication { // for merge into trunk
class SpectralCounts{

 public:

  SpectralCounts();
  ~SpectralCounts();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();

 private:
  // internally-used types
  /**
   * \typedef PeptideSet
   * \brief Collection of peptide objects (not a meta-peptide)
   */
  typedef std::set<PEPTIDE_T*, bool(*)(PEPTIDE_T*, PEPTIDE_T*)> PeptideSet;
  /**
   * \typedef MetaMapping
   * \brief Mapping of peptideSet to MetaProtein
   * Each entry is a set of peptides mapped to a set of proteins of which
   * all contain the set of peptides
   */
  typedef std::map<PeptideSet, MetaProtein, 
              bool(*)(PeptideSet, PeptideSet) > MetaMapping;
  /**
   * \typedef ProteinToPeptides
   * Mapping of Protein objects to a set of peptides that are part
   * of the protein sequence
   */
  typedef std::map<PROTEIN_T*, PeptideSet , 
              bool(*)(PROTEIN_T*, PROTEIN_T*)> ProteinToPeptides;
  /**
   * \typedef MetaToScore
   * \brief Mapping of MetaProtein to the score assigned to it
   */
  typedef std::map<MetaProtein, FLOAT_T, 
              bool(*)(MetaProtein, MetaProtein)> MetaToScore;
  /**
   * \typedef ProteinToMeta
   * \brief Mapping of Protein to MetaProtein to which it belongs
   */
  typedef std::map<PROTEIN_T*, MetaProtein, 
              bool(*)(PROTEIN_T*, PROTEIN_T*)> ProteinToMetaProtein;
  
  // private functions
  void filter_matches(MATCH_COLLECTION_ITERATOR_T* match_collection_it,
                      std::set<MATCH_T*>& match_set);
  void get_peptide_scores(std::set<MATCH_T*>&  matches, 
                          PeptideToScore& peptideToScore);
  void get_protein_scores(PeptideToScore* peptideToScore,
                          ProteinToScore* proteinToScore);
  void get_protein_to_peptides(PeptideToScore* peptideToScore,
                               ProteinToPeptides* proteinToPeptides);
  void get_protein_to_meta_protein(MetaMapping* metaMapping,
                                   ProteinToMetaProtein* proteinToMetaProtein);
  void get_meta_mapping(ProteinToPeptides* proteinToPeptides,
                        MetaMapping& metaMapping);
  void get_meta_ranks(MetaToScore* metaToScore,
                      MetaToRank* metaToRank);
  void get_meta_scores(MetaMapping* metaMapping,
                       ProteinToScore* proteinToScore,
                       MetaToScore* metaToScore);
  void perform_parsimony_analysis(MetaMapping* metaMapping);
  void normalize_peptide_scores(PeptideToScore* peptideToScore);
  void normalize_protein_scores(ProteinToScore* proteinToScore );
  void make_unique_mapping(PeptideToScore* peptideToScore);
  void getSpectra(std::map<std::pair<int,int>, Spectrum*>& spectra);
  int sum_match_intensity(MATCH_T* match,
                          SpectrumCollection* spectra,
                          FLOAT_T bin_width);

  // comparison function declarations
  static bool compare_peptide_sets(PeptideSet, PeptideSet);
  static bool compare_meta_proteins(MetaProtein, MetaProtein);
  static bool sets_are_equal_size(std::pair<PeptideSet, MetaProtein>,
                                  std::pair<PeptideSet, MetaProtein>);
 
}; // class


#endif
