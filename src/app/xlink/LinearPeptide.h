/**
 * \file LinearPeptide.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Object for Defining a Normal peptide in an xlink search
 *****************************************************************************/
#ifndef LINEARPEPTIDE_H_
#define LINEARPEPTIDE_H_

#include "objects.h"
#include "util/utils.h"

#include "XLinkMatch.h"
#include "XLinkBondMap.h"
#include "XLinkablePeptide.h"

#include <vector>

class LinearPeptide : public XLinkMatch {
 protected:

  static std::vector<LinearPeptide> target_linear_peptides_;
  static std::vector<LinearPeptide> decoy_linear_peptides_;
  static void generateAllLinearPeptides(
    std::vector<LinearPeptide>& linear_peptides,
    Database* database,
    PEPTIDE_MOD_T** peptide_mods,
    int num_peptide_mods,
    bool decoy
    );

  Crux::Peptide* peptide_; ///< Peptide this linear peptide referes to
  char* sequence_; ///< sequence of the peptide
  bool is_decoy_; ///< indicator of whether the peptide is a decoy or not
 public:
  /**
   * Default constructor
   */
  LinearPeptide();
  
  /**
   * Constructor with sequence
   */
  LinearPeptide(
    char* sequence_ ///< sequence string
  );
  
  /**
   * Constructor from a Crux peptide
   */
  LinearPeptide(
    Crux::Peptide* peptide ///< peptide object
  );

  /**
   * Default destructor
   */
  virtual ~LinearPeptide() {/*carp(CARP_INFO, "~LinearPeptide");*/};
  
  /**
   *Add candidates to the XLinkMatchCollection that are linear
   */
  static void addCandidates(
    FLOAT_T min_mass, ///< min mass
    FLOAT_T max_mass, ///< max mass
    bool is_decoy,
    XLinkMatchCollection& candidates ///< Vector of possible candidate -inout
    );

  /**
   * returns the candidate type, either a deadlink or a linear candidate
   */
  virtual XLINKMATCH_TYPE_T getCandidateType();
  
  /**
   * \returns the sequence of the peptide in string format
   */
  virtual std::string getSequenceString();
 
  /**
   * \returns the mass of the peptide
   */
  virtual FLOAT_T calcMass(
    MASS_TYPE_T mass_type ///< MONO or AVERAGE
  );
  
  /**
   *\returns a shuffled version of the peptide
   */
  virtual XLinkMatch* shuffle();
  
  /**
   * predicts the ions for this peptide
   */
  virtual void predictIons(
    IonSeries* ion_series, ///< ion series to place the ions
    int charge ///< charge state of the peptide
    );
  
  /**
   * \returns the ion sequence as a string
   */
  std::string getIonSequence(
    Ion* ion ///< ion object
  );
  
  /**
   * \returns the peptide for this match
   */
  virtual Crux::Peptide* getPeptide( 
    int peptide_idx ///< should always be zero
  );
  
  /**
   * \returns the number of missed cleavages
   */
  virtual int getNumMissedCleavages();
  
  /**
   * \returns whether this peptide is modified by a variable mod
   */
  virtual bool isModified();

  static std::vector<LinearPeptide>::iterator getLinearPeptidesBegin(
								     Database* database,
								     PEPTIDE_MOD_T** peptide_mods,
								     int num_peptide_mods,
								     bool decoy
								     );
  static std::vector<LinearPeptide>::iterator getLinearPeptidesEnd(
								   Database* database,
								   PEPTIDE_MOD_T** peptide_mods,
								   int num_peptide_mods,
								   bool decoy);



};

bool compareLinearPeptideMass(const LinearPeptide& xpep1, const LinearPeptide& xpep2);
bool compareLinearPeptideMassToFLOAT(const LinearPeptide& xpep1, FLOAT_T mass);
bool compareLinearPeptideMassToFLOAT2(
  const FLOAT_T& mass,
  const LinearPeptide& pep1);

#endif

/*
 * Local Variables:
 * mode: c 
 * c-basic-offset: 2
 * End:
 */
