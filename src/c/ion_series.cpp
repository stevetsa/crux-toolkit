/*************************************************************************//**
 * \file ion_series.c
 * AUTHOR: Chris Park
 * CREATE DATE: 21 Sep 2006
 * DESCRIPTION: code to support working with a series of ions
 * REVISION: $Revision: 1.52 $
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "objects.h"
#include "ion_series.h"
#include "ion.h"
#include "utils.h"
#include "crux-utils.h"
#include "parameter.h"
#include "peptide.h"
#include "mass.h"
#include "spectrum.h"



/**
 * \struct loss_limit
 * \brief An object that specifies the max amount of neutral loss
 * possible at a given cleavage index. 
 * All numbers are for forward ions(A,B,C) subtract from total to get
 * reverse limit.
 */
struct loss_limit{
  int nh3; ///< the limit to how many NH3 may be lost
  int h2o; ///< the limit to how many H2O may be lost
  // add more if needed for other neutral loss
  // If change this struct, must also modify update_ion_series method
};


/********************************************************
 *ION_SERIES_T
 ********************************************************/
void ION_SERIES_T::init() {
  is_predicted = FALSE;
  peptide = NULL;
  modified_aa_seq = NULL;
  constraint = NULL;
  loss_limit = NULL;
}

/**
 * \returns An (empty) ion_series object.
 */
ION_SERIES_T::ION_SERIES_T() {
  init();
}

/**
 * \brief Creates an ion series for a specific peptide, charge state,
 * and constraint without acutally predicting the ions.
 *
 * Creates a copy of the peptide sequence and calcualtes the peptide
 * mass.
 * Use this method to create ion_series only when few are needed
 * because the memory allocation process is expensive.  Alternatively,
 * use "new_ion_series_generic" and "update_ion_series" in combination
 * to reuse an ion_seires object.
 * \returns A newly allocated ion_series object from the given
 * peptide, charge, and constraint.
 */
ION_SERIES_T::ION_SERIES_T(
  char* peptide, ///< The peptide for this ion series. -in
  int charge, ///< The charge for this ion series -in
  ION_CONSTRAINT_T* constraint ///< constraints which these ions obey.
  )
{
  //carp(CARP_ERROR,"calling init");
  init();
  
  //carp(CARP_ERROR,"copy string");
  // copy the peptide sequence
  this -> peptide = my_copy_string(peptide);
  //carp(CARP_ERROR,"convert_to_mod_aa_seq");
  this -> modified_aa_seq = convert_to_mod_aa_seq(peptide);
  //carp(CARP_ERROR,"calc_sequence_mass");
  this -> peptide_mass = calc_sequence_mass(peptide, MONO);
  this -> charge = charge;
  this -> constraint = constraint;
  this -> peptide_length = strlen(peptide);
  
  // create the loss limit array
  this->loss_limit = 
    (LOSS_LIMIT_T*)mycalloc(this -> peptide_length, sizeof(LOSS_LIMIT_T));
}

/**
 * \brief Creates an ion_series object without a specific peptide.
 * Peptide details are added with the "update_ion_series" method.
 * \returns A newly allocated ion_series object that must be updated
 * for each peptide instance.
 */
ION_SERIES_T::ION_SERIES_T(
  ION_CONSTRAINT_T* constraint, ///< constraints which these ions obey.
  int charge ///< The charge for this ion series -in
  )
{
  //carp(CARP_ERROR,"init");
  init();
  //carp(CARP_ERROR,"constraint");
  this -> constraint = constraint;
  //carp(CARP_ERROR,"charge");
  this -> charge = charge;
  // use max peptide len so loss_limit array can be used for any peptide
  //carp(CARP_ERROR,"loss_limit");
  this -> loss_limit = 
    (LOSS_LIMIT_T*)mycalloc(get_int_parameter("max-length"), 
                            sizeof(LOSS_LIMIT_T));
}

ION_SERIES_T::~ION_SERIES_T() {
  //carp(CARP_ERROR,"ION series destructor");
  //carp(CARP_ERROR,"free(peptide)");
  if(peptide){
    free(peptide);
  }
  //carp(CARP_ERROR,"free modifified_aa_seq");
  if(modified_aa_seq){
    free(modified_aa_seq);
  }

  //carp(CARP_ERROR,"free loss_limit");
  if(loss_limit){
    free(loss_limit);
  }
  // free constraint?

  // iterate over all ions, and free them
  //carp(CARP_ERROR,"free ions");
  for (unsigned int i=0;i<ions.size();i++)
    delete ions[i];
  ions.clear();

  //carp(CARP_ERROR,"done ionseries destructor");
}


/**
 * \brief Updates an ion_series to a specific instance of a peptide
 * sequence. If the ion_series has already generated ions, they will
 * be free-ed. A copy of the peptide sequence is made and all other
 * variables (except ion_constraint) are also updated for the new
 * peptide sequence. 
 */
void ION_SERIES_T::update_ion_series(
  char* peptide, ///< The peptide sequence with no mod characters. -in
  MODIFIED_AA_T* mod_seq ///< modified version of char* sequence -in
  ) 
{
  int ion_type_idx = 0;

  // Initialize the ion_series object for the new peptide sequence
  
  //carp(CARP_ERROR,"free peptide");
  // free old peptide sequence
  if(this->peptide){
    free(this->peptide);
  }
  //carp(CARP_ERROR,"free modified_aa_seq");
  if(this->modified_aa_seq){
    free(this->modified_aa_seq);
  }
  
  //carp(CARP_ERROR,"free ions");
  // iterate over all ions, and free them

  for (unsigned int i=0;i<ions.size();i++)
    delete ions[i];
  ions.clear();
  
  // initialize all num_specific_ion count back to 0
  for(ion_type_idx=0; ion_type_idx < MAX_NUM_ION_TYPE; ++ion_type_idx){
    this->specific_ions[ion_type_idx].clear();
  }
  
  this->is_predicted = FALSE;
  
  // set ion_series for new instance of peptide
  
  // copy the peptide sequence
  this->peptide = my_copy_string(peptide);
  this->peptide_length = strlen(peptide);
  this->modified_aa_seq = copy_mod_aa_seq(mod_seq, 
                                                this->peptide_length);
  
  // Initialize the loss limit array for the new peptide
  int ion_idx;
  for(ion_idx =0; ion_idx < this->peptide_length; ++ion_idx){
    this->loss_limit->h2o = 0;
    this->loss_limit->nh3 = 0;
    // add more initialize count if more added
  }
}

/**
 * Prints a ion_series object to file.
 */
void ION_SERIES_T::print_ion_series(
  FILE* file ///< file for output -out
  )
{
  // check if the ions has been already predicted
  if(!this->is_predicted){
    carp(CARP_ERROR, "ion series has not predicted ions");
    return;
  }
  
  // print header
  fprintf(file, "m/z\tmass\tcharge\tion-series\tpeptide-bond-index\tNH3\tH2O\tISOTOPE\tFLANK\n");
  
  
  // print each ion in the ion series
  int ion_idx;

  
  for(ion_idx=0; ion_idx < numIons(); ++ion_idx){
    ions[ion_idx] -> print_ion(file);
  }
}

/**
 * Prints a ion_series object to file, in GMTK single-ion format.
 */
void ION_SERIES_T::print_ion_series_single_gmtk(
  ION_CONSTRAINT_T* ion_constraint, ///< ion_constraint to obey -in 
  FILE* file,                       ///< file output
  int sentence_idx){

  // create the filtered iterator that will select among the ions
  
  ION_FILTERED_ITERATOR_T ion_iterator;
  // foreach ion in ion iterator, add matched observed peak intensity
  ION_T* ion;
  int frame_idx = 0;

  for (ion_iterator = begin(ion_constraint);
       ion_iterator != end();
       ++ion_iterator) {
    ion = *ion_iterator;
#ifdef BINARY_GMTK
    ion -> print_ion_gmtk_single_binary(file, sentence_idx, frame_idx);
#else
    ion -> print_ion_gmtk_single(file);
    sentence_idx++; // hack to avoid error for not using sentence_idx
#endif
    frame_idx++;
  }
  
  // print a null ion if there are none in this ion series
#ifdef PRINT_NULL_IONS
  for (; frame_idx < MIN_FRAMES; frame_idx++){
#ifdef BINARY_GMTK
    ION_T::print_null_ion_gmtk_single_binary(file, sentence_idx, frame_idx);
#else
    print_null_ion_gmtk_single(file);
    sentence_idx++; // hack to avoid error for not using sentence_idx
#endif
  }
#endif 
}

/**
 * Prints a ion_series object to file, in GMTK paired-ion format.
 */
void ION_SERIES_T::print_ion_series_paired_gmtk(
  ION_CONSTRAINT_T* first_ion_constraint, ///< ion_constraint to obey -in 
  ION_CONSTRAINT_T* second_ion_constraint, ///< ion_constraint to obey -in 
  FILE* file, ///< file output
  int sentence_idx
  )
{
  
  // create the filtered iterator that will select among the ions

  ION_FILTERED_ITERATOR_T ion_iterator;
  // foreach ion in ion iterator, add matched observed peak intensity
  int frame_idx = 0;

  for (ion_iterator = begin(first_ion_constraint);
       ion_iterator != end();
       ++ion_iterator) {
    ION_T* first_ion = *ion_iterator;
    int cleavage_idx = first_ion -> get_ion_cleavage_idx();
    ION_T* second_ion = this -> get_ion_series_ion(second_ion_constraint, cleavage_idx);
    if ( (first_ion == NULL) || (second_ion == NULL) ){
      continue;
    }
    ION_T::print_ion_gmtk_paired_binary(
                                 first_ion, 
                                 second_ion, 
                                 file,
                                 sentence_idx,
                                 frame_idx++);
  }
  
#ifdef PRINT_NULL_IONS
  for (; frame_idx < MIN_FRAMES; frame_idx++){
    ION_T::print_null_ion_gmtk_paired_binary(file, sentence_idx, frame_idx);
  }
#endif
}


/**
 * \brief Find instances of amino acid which can incur neutral
 * losses: H2O (S|T|E|D), NH3(R|K|Q|N).  
 * Set the count of those observed so far for each cleavage index.
 * If no instance of amino acid, the count is assigned to 0
 * The information is used to determine if how many nh3 or h2o neutral
 * losses are possible. 
 */
void ION_SERIES_T::scan_for_aa_for_neutral_loss()
{
  char* sequence = peptide;

  int h2o_aa = 0;
  int nh3_aa = 0;
  LOSS_LIMIT_T* loss_limit_count = NULL; // debug

  // search for the first instance of the amino acids
  int cleavage_idx;
  for(cleavage_idx=0; cleavage_idx < peptide_length; ++cleavage_idx){
    // is the AA  (S|T|E|D) ?
    if(sequence[cleavage_idx] == 'S' ||
       sequence[cleavage_idx] == 'T' ||
       sequence[cleavage_idx] == 'E' ||
       sequence[cleavage_idx] == 'D' )
      {
	//carp(CARP_ERROR,"count");
        loss_limit_count = &this->loss_limit[cleavage_idx];
        loss_limit_count->h2o = ++h2o_aa;
        loss_limit_count->nh3 = nh3_aa;
      }
    // is the AA  (R|K|Q|N) ?
    else if(sequence[cleavage_idx] == 'R' ||
            sequence[cleavage_idx] == 'K' ||
            sequence[cleavage_idx] == 'Q' ||
            sequence[cleavage_idx] == 'N' )
      {
	//carp(CARP_ERROR,"count2");
        loss_limit_count = &loss_limit[cleavage_idx];
        loss_limit_count->nh3 = ++nh3_aa;
        loss_limit_count->h2o = h2o_aa;
      }
    else{
      //carp(CARP_ERROR,"count3");
      loss_limit_count = &loss_limit[cleavage_idx];
      loss_limit_count->h2o = h2o_aa;
      loss_limit_count->nh3 = nh3_aa;
    }
  }
}

/**
 * \brief Creates an array in which element i is the sum of the masses
 * of amino acids 0 to (i-1).  At i=0 is stored the length of the
 * peptide.  
 * \returns an array of ion masses for all sub sequences
 */
FLOAT_T* ION_SERIES_T::create_ion_mass_matrix(
  //char* peptide, ///< The peptide for this ion series. -in
  MODIFIED_AA_T* modified_seq, ///< the sequence
  MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
  int peptide_length ///< the length of the peptide
  )
{
  if( modified_seq == NULL ){
  //if( peptide == NULL ){
    carp(CARP_ERROR, "Cannot create mass matrix from NULL seqence");
    return NULL;
  }
  FLOAT_T* mass_matrix = (FLOAT_T*)mymalloc(sizeof(FLOAT_T)*(peptide_length+1));
  
  // at index 0, the length of the peptide is stored
  mass_matrix[0] = peptide_length;

  // add up AA masses
  int ion_idx = 1;
  // initialize first to be mass of c-term amino acid
  // mass_matrix[ion_idx] = get_mass_amino_acid(peptide[ion_idx-1], mass_type);
  mass_matrix[ion_idx] = get_mass_mod_amino_acid(modified_seq[ion_idx-1], mass_type);
  //++ion_idx;
  //for(; ion_idx <= peptide_length; ++ion_idx){
  for(ion_idx = 2; ion_idx <= peptide_length; ++ion_idx){
    mass_matrix[ion_idx] = mass_matrix[ion_idx-1] + 
      get_mass_mod_amino_acid(modified_seq[ion_idx-1], mass_type);
      //get_mass_amino_acid(peptide[ion_idx-1], mass_type);
  }
  // DEBUGGING
  /*
  fprintf(stderr, "seq:");
  for(ion_idx = 0; ion_idx < peptide_length; ++ion_idx){
    fprintf(stderr, "\t%s", modified_aa_to_string(modified_seq[ion_idx]));
  }
  fprintf(stderr, "\nmas:");
  for(ion_idx = 0; ion_idx < peptide_length; ++ion_idx){
    fprintf(stderr, "\t%.2f", get_mass_mod_amino_acid(modified_seq[ion_idx], MONO));
  }
  fprintf(stderr, "\nsum:");
  for(ion_idx = 0; ion_idx < peptide_length; ++ion_idx){
    fprintf(stderr, "\t%.2f", mass_matrix[ion_idx+1]);
  }
  fprintf(stderr, "\n");
  */
  return mass_matrix;
}

/**
 * user must ensure that there is enough space for this ion
 * adds ion to ion_series' master ion_array and if B|Y ion to the specific ion_array
 */
void ION_SERIES_T::add_ion_to_ion_series(
  ION_T* ion ///< ion to add -in
  )
{
  // add ion to ion series
  ions.push_back(ion);
  
  // add a pointer of ion to the specific ion_type array

  specific_ions[ion -> get_ion_type()].push_back(ion);
}

/**
 * helper function: add_ions
 * add all the ions to ion_series up to the max charge
 *\returns TRUE if successfully adds all ions, else FALSE
 */
BOOLEAN_T ION_SERIES_T::add_ions_by_charge(
  FLOAT_T mass, ///< the base mass of the ion to add
  int cleavage_idx, ///< the absolute cleavage index (A,B,C from left X,Y,Z from right)
  ION_TYPE_T ion_type ///< the ion type of the ions to be added
  )
{
  int charge_idx = 1;
  ION_T* ion = NULL;
  int max_charge;

  // check if there's enough space to add the new ions
  if(numIons() + constraint->get_max_charge() > MAX_IONS){
    carp(CARP_FATAL, "exceeds ion array size in ion_series");
    return FALSE;
  }
  
  // set the max charge, the maximum cannot exceed the precursor ion's charge
  if(constraint->get_max_charge() > charge){
    max_charge = this->charge;
  }
  else{
    max_charge = constraint -> get_max_charge();
  }

  // iterate over all different charge
  for(; charge_idx <= max_charge; ++charge_idx){
    // create ion
    ion = new ION_T(ion_type, 
		    cleavage_idx, 
		    charge_idx, 
		    this->peptide, 
		    constraint->get_mass_type(), 
		    mass); 
    // add ion to ion series
    this -> add_ion_to_ion_series(ion);
  }
  
  return TRUE;
}


/**
 * Creates all the ions with no modifications up to the max charge
 * Adds each ion to ion_series
 *\returns TRUE if successfully generates all the ions, else FALSE
 */
BOOLEAN_T ION_SERIES_T::generate_ions_no_modification(
  FLOAT_T* mass_matrix ///< the mass matrix that stores the mass
  )
{
  if( mass_matrix == NULL ){
    carp(CARP_ERROR,
         "Cannot generate ions from NULL ion series or mass matrix");
    return FALSE;
  }
  int cleavage_idx = 1;
  FLOAT_T mass = 0;

  // get peptide length
  int peptide_length = (int)mass_matrix[0];

  // iterate over all cleavage index
  for(; cleavage_idx < peptide_length; ++cleavage_idx){
    
    // add A ion
    if(constraint-> get_ion_type() == A_ION 
       || constraint-> get_ion_type() == BYA_ION 
       || constraint-> get_ion_type() == ALL_ION){

      // set mass
      mass = mass_matrix[cleavage_idx];
      
      if(constraint->get_mass_type() == MONO){
        mass -= MASS_CO_MONO; 
      }
      else{ // average
        mass -= MASS_CO_AVERAGE; 
      }
      
      // add ions up to max charge
      if(!this -> add_ions_by_charge(mass, cleavage_idx, A_ION)){
        carp(CARP_ERROR, "failed to add ions by different charge for A ion");
      return FALSE;
      }
    }
    
    // add B ion
    if(constraint -> get_ion_type() == ALL_ION 
       || constraint -> get_ion_type() == BY_ION
       || constraint-> get_ion_type() == BYA_ION
       || constraint-> get_ion_type() == B_ION){
      
      // set mass
      mass = mass_matrix[cleavage_idx];
      
      // add ions up to max charge
      if(!this -> add_ions_by_charge(mass, cleavage_idx, B_ION)){
        carp(CARP_ERROR, "failed to add ions by different charge for B ion");
        return FALSE;
      }
    }
    
    // add C ion
    if(constraint -> get_ion_type() == C_ION || constraint -> get_ion_type() == ALL_ION){
      // set mass
      mass = mass_matrix[cleavage_idx];
      
      if(constraint->get_mass_type() == MONO){
        mass += MASS_NH3_MONO; 
      }
      else{ // average
        mass += MASS_NH3_AVERAGE; 
      }
      
      // add ions up to max charge
      if(!this -> add_ions_by_charge(mass, cleavage_idx, C_ION)){
        carp(CARP_ERROR, "failed to add ions by different charge for C ion");
        return FALSE;
      }
    }
    
    // add X ion
    if(constraint -> get_ion_type() == X_ION || constraint-> get_ion_type() == ALL_ION){
      // set mass 
      mass = mass_matrix[(int)mass_matrix[0]] - mass_matrix[(int)mass_matrix[0] - cleavage_idx];

      if(constraint -> get_mass_type() == MONO){
        mass += MASS_CO_MONO + MASS_H2O_MONO;     
      }
      else{ // average
        mass += MASS_CO_AVERAGE + MASS_H2O_AVERAGE; 
      }
      
      // add ions up to max charge
      if(!this -> add_ions_by_charge(mass, cleavage_idx, X_ION)){
        carp(CARP_ERROR, "failed to add ions by different charge for X ion");
        return FALSE;
      }
    }
    
    // add Y ion
    if(constraint->get_ion_type() == ALL_ION || 
       constraint->get_ion_type() == BY_ION ||
       constraint->get_ion_type() == BYA_ION ||
       constraint->get_ion_type() == Y_ION){

      // set mass 
      mass = mass_matrix[(int)mass_matrix[0]] - mass_matrix[(int)mass_matrix[0] - cleavage_idx];
      
      if(constraint->get_mass_type() == MONO){
        mass += MASS_H2O_MONO; 
      }
      else{ // average
        mass += MASS_H2O_AVERAGE;
      }
      
      // add ions up to max charge
      if(!this -> add_ions_by_charge(mass, cleavage_idx, Y_ION)){
        carp(CARP_ERROR, "failed to add ions by different charge Y ion");
        return FALSE;
      }
      
    }
    
    // add Z ion
    if(constraint -> get_ion_type() == Z_ION ||
       constraint -> get_ion_type() == ALL_ION ){

      // set mass 
      mass = mass_matrix[(int)mass_matrix[0]] - mass_matrix[(int)mass_matrix[0] - cleavage_idx];
      
      if(constraint -> get_mass_type() == MONO){
        mass = mass - MASS_NH3_MONO + MASS_H2O_MONO; 
      }
      else{ // average
        mass = mass - MASS_NH3_AVERAGE + MASS_H2O_AVERAGE;
      }
      
      // add ions up to max charge
      if(!this -> add_ions_by_charge(mass, cleavage_idx, Z_ION)){
        carp(CARP_ERROR, "failed to add ions by different charge Z ion");
        return FALSE;
      }
    }
  }

  // add P ion(precursor ion)?
  if(constraint->get_precursor_ion()){
    
    // set mass 
    mass = mass_matrix[(int)mass_matrix[0]];
    
    // mass type
    if(constraint->get_mass_type() == MONO){
      mass += MASS_H2O_MONO; 
    }
    else{ // average
      mass += MASS_H2O_AVERAGE;
    }
    
    // add ions up to max charge
    if(!add_ions_by_charge(mass, (int)mass_matrix[0], P_ION)){
      carp(CARP_ERROR, "failed to add ions by different charge P ion");
      return FALSE;
    }
  }
  return TRUE;
}

/**
 * The modification depends on the loss/add && if the ion contains RKQN or STED
 * The number of losses possible cannot exceed the number of RKQN or STED in the ion
 * The loss_limit array in the ion_series must be populated prior to this method call
 *\returns TRUE if the ion can lose the mod_type modification, else FALSE
 */
BOOLEAN_T ION_SERIES_T::can_ion_lose_modification(
  ION_T* ion, ///< the ion to check if can lose nh3 -in
  ION_MODIFICATION_T mod_type, ///< generate ions of this modification_type -in/out
  int increment  ///< the add/loss of the modification
  )
{
  int cleavage_idx = ion -> get_ion_cleavage_idx();

  // check for NH3 modification
  if(mod_type == NH3){
    // adding is ok
    if(increment >= 0){
      return TRUE;
    }
    
    // is forward ion_type(ABC)?
    if(ion -> is_forward_ion_type()){
      // does this ion contain enough counts of the RKQN
      if(-increment >  (&loss_limit[cleavage_idx-1])->nh3){
        return FALSE;
      }
      return TRUE;
    }
    else{// backward ions XYZ
      // does this ion contain enough counts of the RKQN
      if(cleavage_idx == peptide_length){
        if(-increment > (&loss_limit[this->peptide_length-1])->nh3){
          return FALSE;
        }
      }
      else if(-increment >  
              ((&loss_limit[peptide_length-1])->nh3 - 
               (&loss_limit[peptide_length - cleavage_idx - 1])->nh3)){
        return FALSE;
      }
      return TRUE;
    }
  }

  // check for H2O modification
  if(mod_type == H2O){
    // adding is ok
    if(increment >= 0){
      return TRUE;
    }
    
    // is forward ion_type(ABC)?
    if(ion -> is_forward_ion_type()){
      // does this ion contain enough counts of the STED
      if(-increment >  (&loss_limit[cleavage_idx-1])->h2o){
        return FALSE;
      }
      return TRUE;
    }
    else{// backward ions XYZ
      // does this ion contain enough counts of the STED
      if(cleavage_idx == this->peptide_length){
        if(-increment > (&loss_limit[this->peptide_length-1])->h2o){
          return FALSE;
        }
      }
      else if(-increment >  
              ((&loss_limit[peptide_length-1])->h2o - 
               (&loss_limit[peptide_length - cleavage_idx - 1])->h2o)){
        return FALSE;
      }
      return TRUE;
    }
  }
  
  // check for ISOTOPE modification
  else if(mod_type == ISOTOPE){
    // adding is ok
    if(increment >= 0 && ion -> get_ion_type() != P_ION){
      return TRUE;
    }
    // add more constraint if needed
  }
  // check for FLANK modification
  else if(mod_type == FLANK){
    // only add flanking ions to type B,Y ions
    if(ion -> get_ion_type() == B_ION || ion -> get_ion_type() == Y_ION){
      // only add ions with no modifications
      if(!ion -> ion_is_modified()){
        return TRUE;
      }
    }
  }
  
  return FALSE;
}

/**
 * creates all the ions with specific modifications up to the max charge
 * copies all the existing ions that can be modified,
 * then applies the different modifications then adds the new modified ions to ion_series
 *\returns TRUE if successfully generates all the ions with modifications, else FALSE
 */
BOOLEAN_T ION_SERIES_T::generate_ions(
  ION_MODIFICATION_T mod_type ///< generate ions of this modification_type -in/out
  )
{
  int ion_idx = 0;
  int total_ion = numIons();
  ION_T* working_ion = NULL;
  ION_T* new_ion = NULL;
  int* modifications = constraint->get_modifications();

  // modification index
  int type_idx = 0;
  int type_increment = 1;

  // check if there's enough space to add the new ions
  if(numIons()*2 > MAX_IONS){
    carp(CARP_FATAL, "exceeds ion array size in ion_series");
    return FALSE;
  }

  // reset modification increment, if mod_type loss
  if(modifications[mod_type] < 0){
    type_increment = -1;
  }
  
  // iterate over all the ions to determine which ones should be copies and modified
  for(; ion_idx < total_ion; ++ion_idx){
    working_ion = ions[ion_idx];
    
    // can this ion generate a mod_type modification?, if not skip to next ion
    if(!(can_ion_lose_modification(working_ion, mod_type, type_increment))){      
      continue;
    }
     
    // add/sub thorugh all mod_type modifications!!!
    for(type_idx = type_increment; abs(type_idx) <= abs(modifications[mod_type]); ){
      // copy the src ion, into new ion
      new_ion = new ION_T();
      ION_T::copy_ion(working_ion, new_ion, working_ion -> get_ion_peptide_sequence());
      
      // add the modification to the new ion
      new_ion -> add_modification(mod_type, type_idx, constraint -> get_mass_type());
      
      // add ion to ion_series
      add_ion_to_ion_series(new_ion);
     
      // can this ion generate a mod_type modification for the next count of modification?, 
      if(!(can_ion_lose_modification(working_ion, mod_type, 
              (type_idx += type_increment)))){
        break;
      }
    }
  }
  return TRUE;
}

/**
 * creates all the flanking ions up to the max charge
 * can only create flanking ions that are B|Y ions and don't have modification
 * assumes the ions with no modification all are at the begining of the ion[] in ion_series
 * copies all the existing ions that can be modified,
 * then applies the different modifications then adds the new modified ions to ion_series
 *\returns TRUE if successfully generates all the ions with modifications, else FALSE
 */
BOOLEAN_T ION_SERIES_T::generate_ions_flank()
{
  int ion_idx = 0;
  int total_ion = numIons();
  ION_T* working_ion = NULL;
  ION_T* new_ion = NULL;
  ION_T* new_ion2 = NULL;
  int* modifications = constraint -> get_modifications();

  // modification index
  int type_idx = 0;
  int type_increment = 1;

  // check if there's enough space to add the new ions
  if(numIons()*2 > MAX_IONS){
    carp(CARP_FATAL, "exceeds ion array size in ion_series");
    return FALSE;
  }
  
  // iterate over all the ions to determine which ones should be copies and modified
  for(; ion_idx < total_ion; ++ion_idx){
    working_ion = ions[ion_idx];
    
    // no more ions that are not modified, thus done
    if(working_ion -> get_ion_single_modification_count(NH3) != 0){
      break;
    }

    // can this ion generate a mod_type modification?, if not skip to next ion
    if(!can_ion_lose_modification(working_ion, FLANK, type_increment)){      
      continue;
    }
     
    // add/sub thorugh all mod_type modifications!!!
    for(type_idx = type_increment; type_idx <= modifications[FLANK]; type_idx += type_increment){
      // copy the src ion, into new ion
      new_ion = new ION_T();
      new_ion2 = new ION_T();
      ION_T::copy_ion(working_ion, new_ion, working_ion -> get_ion_peptide_sequence());
      ION_T::copy_ion(working_ion, new_ion2, working_ion -> get_ion_peptide_sequence());
      
      // add the modification to the new ion
      new_ion -> add_modification(FLANK, type_idx, constraint -> get_mass_type());
      new_ion2 -> add_modification(FLANK, -type_idx, constraint -> get_mass_type());

      // add ion to ion_series
      add_ion_to_ion_series(new_ion);
      add_ion_to_ion_series(new_ion2);
    }
  }
  return TRUE;
}

/**
 * \brief The engine of ion series. Predicts all the ions from the
 * peptide that meet the ion constraint. All predicted ions are stored
 * in the ion_series as ion objects. 
 */
void ION_SERIES_T::predict_ions()
{
  if(is_predicted){
    carp(CARP_WARNING, "The ion series has already been predicted and added");
    return;
  }

  //carp(CARP_ERROR,"Predicting mass matrix");
  // create a mass matrix
  FLOAT_T* mass_matrix = 
    create_ion_mass_matrix(modified_aa_seq, constraint -> get_mass_type(), peptide_length);  
  /*
  printf("cumulative mass sum is:\n");
  int idx = 0;
  for(idx = 0; idx < mass_matrix[0]; idx++){
    printf("%i\t%f\n", idx, mass_matrix[idx]);
  }
  */
  // scan for the first and last  (S, T, E, D) and (R, K, Q, N), 
  // initialize to determine modification is ok.
  // the first, last of STED, RKQN are stored in ion_series.

  //carp(CARP_ERROR,"Scanning for neutral losses");
  scan_for_aa_for_neutral_loss();
  
  //carp(CARP_ERROR,"generating ions no modification");
  // generate ions without any modifications
  if(!generate_ions_no_modification(mass_matrix)){
    carp(CARP_FATAL, "failed to generate ions, no modifications");
  }

  // create modification ions?
  if(constraint-> get_use_neutral_losses()){
    
    // generate ions with nh3 modification
    if(abs(constraint-> get_modification(NH3)) > 0){
      if(!generate_ions(NH3)){
        carp(CARP_FATAL, "failed to generate ions, NH3 modifications");
      }
    }
    
    // generate ions with h2o modification
    if(abs(constraint-> get_modification(H2O)) > 0){
      if(!generate_ions(H2O)){
        carp(CARP_FATAL, "failed to generate ions, H2O modifications");
      }
    }

    // generate ions with isotope modification
    if(constraint -> get_modification(ISOTOPE) > 0){
      if(!generate_ions(ISOTOPE)){
        carp(CARP_FATAL, "failed to generate ions, ISOTOPE modifications");
      }
    }

    // generate ions with flank modification
    if(constraint-> get_modification(FLANK) > 0){
      if(!generate_ions_flank()){
        carp(CARP_FATAL, "failed to generate ions, FLANK modifications");
      }
    }
    
    // add more modifications here

  }
  
  // ion series now been predicted
  is_predicted = TRUE;
  
  // free mass matrix
  free(mass_matrix);
}

/**
 * Assign peaks to the nearest ions, within a tolerance (set in param file)
 */
void ION_SERIES_T::ion_series_assign_nearest_peaks(
    SPECTRUM_T* spectrum){

  //  FLOAT_T max = 0.5; // TODO set in param file 
  FLOAT_T max = get_double_parameter("ion-tolerance"); 
  ION_T* ion = NULL;

  ION_ITERATOR_T iterator;
  PEAK_T* peak = NULL;

  for(iterator = ions.begin();
      iterator != ions.end();
      ++iterator) {
    ion = *iterator;
    FLOAT_T mz = ion -> get_ion_mass_z(); // TODO change to mz, not mass_z
    peak = spectrum -> get_nearest_peak(mz, max);
    ion -> set_ion_peak(peak);
  }
}

/**
 * Copies ion_series object from src to dest.
 *  must pass in a memory allocated ION_SERIES_T* dest
 * does not copy the loss_limit.
 */
void ION_SERIES_T::copy_ion_series(
  ION_SERIES_T* src,///< ion to copy from -in
  ION_SERIES_T* dest///< ion to copy to -out
  )
{
  ION_T* src_ion = NULL;
  ION_T* dest_ion = NULL;
  
  dest->peptide = my_copy_string(src->peptide);
  dest->charge = src->charge;
  dest->peptide_length = src->peptide_length;
  //mod seq???

  // add copy of pointer ion constraint
  dest->constraint = src->constraint;

  // add copy ion, add ion_filtered_iterator

  ION_ITERATOR_T iterator;
  // iterate over all ions in src and copy them into dest
  for (iterator = src -> ions.begin();
       iterator != src -> ions.end();
       ++iterator) {
    src_ion = *iterator;
    // add ion
    dest_ion = new ION_T();
    ION_T::copy_ion(src_ion, dest_ion, dest->peptide);
    dest -> add_ion_to_ion_series(dest_ion);
  }

  dest->is_predicted = TRUE;
}

/*************************************
 * ION_SERIES_T: get and set methods
 ************************************/

/**
 * \returns the ion that meets the criteria or NULL
 * TODO possibly reimplement if linear scan is too slow
 */
ION_T* ION_SERIES_T::get_ion_series_ion(
    ION_CONSTRAINT_T* ion_constraint,
    int cleavage_idx
    ){

  ION_FILTERED_ITERATOR_T ion_iterator;
  
  ION_T* ion = NULL;

  for (ion_iterator = begin(ion_constraint);
       ion_iterator != end();
       ++ion_iterator) {
    ion = *ion_iterator;
    if(ion -> get_ion_cleavage_idx() == cleavage_idx){
      return ion;
    }
  }
  return NULL;
}

/**
 *\returns the peptide length of which the ions are made
 */
int ION_SERIES_T::get_ion_series_peptide_length()
{
  return peptide_length;
}

/**
 * User should not free the peptide sequence seperate from the ion_series
 *\returns a pointer to the original parent peptide sequence of the ion_series object
 */
char* ION_SERIES_T::get_ion_series_peptide()
{
  return peptide;
}

/**
 * copies in the peptide sequence to heap allocated sequence.
 * set the parent peptide sequence of the ion_series object
 */
void ION_SERIES_T::set_ion_series_peptide(
  char* peptide///< the peptide sequence to set -in
  )
{
  // free previous sequence
  if(this -> peptide != NULL){
    free(this -> peptide);
  }
  this -> peptide = my_copy_string(peptide);
}

/**
 *\returns the charge of the ion_series object
 */
int ION_SERIES_T::get_ion_series_charge()
{
  return charge;
}

/**
 * set the charge of the ion_series object
 */
void ION_SERIES_T::set_ion_series_charge(
  int charge///< the charge of the ion -in
  )
{
  this->charge = charge;
}

/**
 *\returns the constraint of the ion_series object
 */
ION_CONSTRAINT_T* ION_SERIES_T::get_ion_series_ion_constraint()
{
  return constraint;
}

/**
 * set the of the ion_series object
 */
void ION_SERIES_T::set_ion_series_ion_constraint(
  ION_CONSTRAINT_T* constraint///<  -in
  )
{
  this->constraint = constraint;
}

/**
 *\returns the total number of ions in the ion_series object
 */
int ION_SERIES_T::numIons()
{
  return ions.size();
}

/**
 *\returns the total number of ion_type in the ion_series object
 */
int ION_SERIES_T::numSpecificIons(
  ION_TYPE_T ion_type ///< the type of ions -in
  )
{
  return specific_ions[ion_type].size();
}

ION_ITERATOR_T ION_SERIES_T::begin() {
  return ions.begin();
}

ION_ITERATOR_T ION_SERIES_T::end() {
  return ions.end();
}



/*************************
 * ION_CONSTRAINT methods
 *************************/


/**
 *\returns an empty heap allocated ion_constraint
 */
void ION_CONSTRAINT_T::init(){
  int modification_idx = 0;
  // initialize all modifications count to 0
  for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
    modifications[modification_idx] = 0;
  }
  pointer_count = 1;
  min_charge = 0;
  use_neutral_losses = FALSE;
  exact_modifications = FALSE;
  
}

void ION_CONSTRAINT_T::init(
  MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
  int max_charge, ///< max charge of the ions <= the parent peptide's charge
  ION_TYPE_T ion_type, ///< the ion types the peptide series should include
  BOOLEAN_T precursor_ion  ///< should include precursor ion?
  ) 
{
  //set all fields of constraint
  init();
  this -> mass_type = mass_type;
  this -> max_charge = max_charge;
  this -> ion_type = ion_type;
  this -> precursor_ion = precursor_ion;
}

ION_CONSTRAINT_T::ION_CONSTRAINT_T() {
  init();
}

/**
 * modification, all modifications 0
 * add more modifications as needed using the set_ion_constraint_modification
 *\returns a new heap allocated ion_constraint
 */
ION_CONSTRAINT_T::ION_CONSTRAINT_T(
  MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
  int max_charge, ///< max charge of the ions <= the parent peptide's charge
  ION_TYPE_T ion_type, ///< the ion types the peptide series should include
  BOOLEAN_T precursor_ion  ///< should include precursor ion?
  )
{
  init(mass_type, max_charge, ion_type, precursor_ion);
}

/**
 * \brief Create a new ion constraint based on the score type and the
 * charge of the peptide to be modeled.  Uses other
 * new_ion_constraint_ methods for some types.
 *
 * \returns A newly allocated ion constraint.
 */
ION_CONSTRAINT_T::ION_CONSTRAINT_T(
  SCORER_TYPE_T score_type,
  int charge
){
  
  switch(score_type){
  case SP:
    init_ion_constraint_sequest_sp(charge);
    break;
  case XCORR:
    init_ion_constraint_sequest_xcorr(charge);
    break;
  case DOTP:
  case LOGP_EXP_SP:
    //case LOGP_BONF_EXP_SP:
    //case LOGP_EVD_XCORR:
  case DECOY_XCORR_QVALUE:
  case DECOY_PVALUE_QVALUE:
  case LOGP_BONF_EVD_XCORR:
  case LOGP_WEIBULL_SP:
  case LOGP_BONF_WEIBULL_SP:
  case LOGP_WEIBULL_XCORR:
  case LOGP_BONF_WEIBULL_XCORR:
  case Q_VALUE:
  case PERCOLATOR_SCORE:
  case LOGP_QVALUE_WEIBULL_XCORR:
      init(get_mass_type_parameter("fragment-mass"),
	   charge,
	   get_ion_type_parameter("primary-ions"),
	   get_boolean_parameter("precursor-ions")); 
    break;
  }
}

/**
 * modification, sets all fields for gmtk settings
 *\returns a new heap allocated ion_constraint
 */
void ION_CONSTRAINT_T::init_ion_constraint_gmtk(
  int charge 
  )
{
  int max_charge = 1;
  if(charge > 1){
    max_charge = charge;
  }  
  init(MONO, max_charge, ALL_ION, FALSE);

  // set all modifications count for gmtk
  use_neutral_losses = TRUE;
  min_charge = 1;
  modifications[NH3] = -1;
  modifications[H2O] = -1;
  modifications[ISOTOPE] = 0;
  modifications[FLANK] = 0;
}

ION_CONSTRAINT_T* ION_CONSTRAINT_T::new_gmtk(
    int charge ///< the charge of the peptide for which to predict ions
  )
{
  ION_CONSTRAINT_T* ans = new ION_CONSTRAINT_T();
  ans -> init_ion_constraint_gmtk(charge);
  return ans;
}

/**
 * modification, sets all fields for sequest settings
 *\returns a new heap allocated ion_constraint
 */
void ION_CONSTRAINT_T::init_ion_constraint_sequest(
  MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
  int max_charge, ///< the maximum charge of the ions. 
                  ///< cannot exceed the parent peptide's charge
  ION_TYPE_T ion_type, ///< the ion types the peptide series should include
  BOOLEAN_T precursor_ion ///< should include precursor ion?
  )
{
  init(mass_type, max_charge, ion_type, precursor_ion);

  // set                                                     
  use_neutral_losses = TRUE;
  
  // set all modifications count for sequest
  modifications[NH3] = 1;
  modifications[H2O] = 1;
  modifications[ISOTOPE] = 1;
  modifications[FLANK] = 1;
}

/**
 * modification, sets all fields for sequest Sp scoring settings
 * make B, Y type ions
 *\returns a new heap allocated ion_constraint
 */
void ION_CONSTRAINT_T::init_ion_constraint_sequest_sp(
  int charge ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
  )
{
  // charge = 1;
  if(charge == 1){
    init(MONO, 1, BY_ION, FALSE);
  }  
  else{
    --charge;
    init(MONO, charge, BY_ION, FALSE);
  }
  
  // set                                                     
  use_neutral_losses = TRUE;
  
  // set all modifications count for sequest
  modifications[NH3] = 0;
  modifications[H2O] = 0;
  modifications[ISOTOPE] = 0;
  modifications[FLANK] = 0;
}


/**
 * modification, sets all fields for Sequest Xcorr scoring settings
 * make B, Y, A type ions
 *\returns a new heap allocated ion_constraint
 */
void ION_CONSTRAINT_T::init_ion_constraint_sequest_xcorr(
  int charge ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
  )
{
  // charge = 1;
  if(charge == 1){
    init(MONO, 1, BYA_ION, FALSE);
  }  
  else{
    --charge;
    init(MONO, charge, BYA_ION, FALSE);
  }
  
  // set                                                     
  use_neutral_losses = TRUE;
  
  // set all modifications count for sequest
  modifications[NH3] = 0;// -1;
  modifications[H2O] = 0;// -1;
  modifications[ISOTOPE] = 0;// not sure
  modifications[FLANK] = 0;
}


/**
 * Frees an allocated ion_constraint object.
 */
ION_CONSTRAINT_T::~ION_CONSTRAINT_T() {

}

/**
 * Frees an allocated ion_constraint object.
 */
void ION_CONSTRAINT_T::free(ION_CONSTRAINT_T* ion_constraint) 
{
  //carp(CARP_ERROR,"Inside ion_constraint free");
  ion_constraint -> pointer_count--;
  if (ion_constraint -> pointer_count == 0){
    //carp(CARP_ERROR,"really deleting ion_constraint");
    delete ion_constraint;
  }
  //carp(CARP_ERROR,"Done ion_constraint free");
}


/**
 * copies ion_constraint pointer
 */
ION_CONSTRAINT_T* ION_CONSTRAINT_T::copy_ion_constraint_ptr(ION_CONSTRAINT_T* constraint)
{
  constraint -> pointer_count++;
  return constraint;
}


/**
 * copies ion_constraint object from src to dest
 * must pass in a memory allocated ION_CONSTRAINT_T dest
 */
void ION_CONSTRAINT_T::copy_ion_constraint(
  ION_CONSTRAINT_T* src,///< ion_constraint to copy from -in
  ION_CONSTRAINT_T* dest///< ion_constraint to copy to -out
  )
{
  int modification_idx = 0;
  dest->use_neutral_losses = src->use_neutral_losses;

  // if use natural loss, copy
  if(src->use_neutral_losses){
    // iterate over all modifications a update new constraint
    for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
      dest->modifications[modification_idx] = src->modifications[modification_idx];
    }
  }
  
  dest->mass_type = src->mass_type;
  dest->max_charge = src->max_charge;
  dest->ion_type = src->ion_type;
  dest->precursor_ion = src->precursor_ion;
}

/** FIXME!!!! double check
 * Determines if a ion satisfies a ion_constraint.
 * \returns TRUE if the constraint is satisified. FALSE if not.
 */
BOOLEAN_T ION_CONSTRAINT_T::is_satisfied(
   ION_T* ion ///< query ion -in
   )
{
  int* counts = NULL;

  // TODO Fix
  BOOLEAN_T return_val = TRUE;
  // print_ion(ion, stderr);
  // fprintf(stderr, "%i->%i\n", ion_constraint->min_charge, ion_constraint->max_charge);
  // check ion type
  ION_TYPE_T ion_type = ion -> get_ion_type();
  if(
     !(ion_type == this ->ion_type)
      
     &&
     
     !((this->ion_type == BY_ION) && (ion_type == B_ION || ion_type == Y_ION)) 
     
     &&

     !((this->ion_type == BYA_ION) 
          && 
       (ion_type == B_ION || ion_type == Y_ION || ion_type == A_ION)) 
     
     &&
     
     !(this->ion_type == ALL_ION)

     ){
     
    // precursor ion?
    if(!(this->precursor_ion && ion_type == P_ION)){
      return_val = FALSE;
    }
  }
  
  // check charge
  if(ion -> get_ion_charge() > this->max_charge){
    return_val = FALSE;
  }
  
  if(ion -> get_ion_charge() < this->min_charge){
    return_val = FALSE;
  }

  // check modifications
  counts = ion -> get_ion_modification_counts();
  int mod_idx;
  for(mod_idx=0; mod_idx < MAX_MODIFICATIONS; ++mod_idx){
    if(this->modifications[mod_idx] >= 0){
      if(counts[mod_idx] > this->modifications[mod_idx]){
        return_val = FALSE;
        break;
      }
    }
    else{
      if(counts[mod_idx] < this->modifications[mod_idx]){
        return_val = FALSE;
        break;
      }
    }
    if (this->exact_modifications){
      if(counts[mod_idx] != this->modifications[mod_idx]){
        return_val = FALSE; 
        break;
      }
    }
  }
  
  // Add more checks here as more contraints are added

  // fprintf(stderr, "r = %i\n", return_val);
  return return_val;
}


/**
 * sets the modification count
 * can only add isotopes
 */
void ION_CONSTRAINT_T::set_ion_constraint_modification(
  ION_MODIFICATION_T mod_type, ///< ion modification type -in
  int count  ///< the count of the modification -in  
  )
{
  // set modification count, can only add for isotope
  if(mod_type != ISOTOPE){
    modifications[mod_type] = count;
  }
  else{
    modifications[mod_type] = abs(count);
  }

  // set neutral loss to TRUE if needed
  if(!use_neutral_losses){
    use_neutral_losses = TRUE;
  }
}

/**
 * sets the exact modification boolean 
 */
void ION_CONSTRAINT_T::set_ion_constraint_exactness(
  BOOLEAN_T exactness ///< whether to use exact mods or not -in
  ){
  exact_modifications = exactness;
  if (exactness == TRUE){
    min_charge = max_charge;
  }
}
 
/**
 * gets the modification count for specific mod_type
 */
int ION_CONSTRAINT_T::get_modification(
  ION_MODIFICATION_T mod_type ///< ion modification type -in
  )
{
  return modifications[mod_type];
}

/**
 * gets the modification count array for specific mod_types
 */
int* ION_CONSTRAINT_T::get_modifications() {
  return modifications;
}

/**
 * gets the mass type of the ion_constraint
 */
MASS_TYPE_T ION_CONSTRAINT_T::get_mass_type()
{
  return mass_type;
}

/**
 * get maximum charge of the ions, cannot exceed the parent peptide's charge
 */
int ION_CONSTRAINT_T::get_max_charge()
{
  return max_charge;
}

/**
 * get the ion types the peptide series should include
 */
ION_TYPE_T ION_CONSTRAINT_T::get_ion_type() 
{
  return ion_type;
}

/*
 * get whether a precursor-ion satisfies this constraint
 */
BOOLEAN_T ION_CONSTRAINT_T::get_precursor_ion() 
{
  return precursor_ion;
}

/*
 * get if ions should include neutral losses
 */
BOOLEAN_T ION_CONSTRAINT_T::get_use_neutral_losses()
{
  return use_neutral_losses;
}

/**********************************
 * ION_FILTERED_ITERATOR_T object
 **********************************/
/**
 * Only copies in the constraint as pointer
 * Instantiates a new ion_filtered_iterator object from ion_series.
 * \returns a ION_FILTERED_ITERATOR_T object.
 */

ION_FILTERED_ITERATOR_T ION_SERIES_T::begin(ION_CONSTRAINT_T* constraint) {
  return FilteredIterator(this, constraint);
}

ION_FILTERED_ITERATOR_T::FilteredIterator() {
  ion_series = NULL;
  constraint = NULL;
}

ION_FILTERED_ITERATOR_T::FilteredIterator(ION_SERIES_T* ion_series, ION_CONSTRAINT_T* constraint) {
  this -> ion_series = ion_series;
  this -> constraint = constraint;
  
  if(constraint-> get_ion_type() == ALL_ION ||
     constraint-> get_ion_type() == BY_ION ||
     constraint-> get_ion_type() == BYA_ION){
    current = ion_series -> begin();
    end_iter = ion_series -> end();
  }
  else{
    current = ion_series -> specific_ions[constraint->get_ion_type()].begin();
    end_iter = ion_series -> specific_ions[constraint->get_ion_type()].end();
  }
  satisfyConstraint();
}


ION_FILTERED_ITERATOR_T::~FilteredIterator() {
  ;
}

ION_FILTERED_ITERATOR_T& ION_FILTERED_ITERATOR_T::operator=(const ION_FILTERED_ITERATOR_T& other) {
  constraint = other.constraint;
  current = other.current;
  ion_series = other.ion_series;
  end_iter = other.end_iter;
  return (*this);
}
  
bool ION_FILTERED_ITERATOR_T::operator==(const ION_FILTERED_ITERATOR_T& other) {
  return (constraint == other.constraint) && (current == other.current); 
}

bool ION_FILTERED_ITERATOR_T::operator!=(const ION_FILTERED_ITERATOR_T& other) {
  return (constraint != other.constraint) || (current != other.current);
}

bool ION_FILTERED_ITERATOR_T::operator!=(const ION_ITERATOR_T& other) {
  return (current != other);
}

ION_FILTERED_ITERATOR_T& ION_FILTERED_ITERATOR_T::operator++() {
  if (current != 
      ion_series -> end()) { 
      ++current;
      satisfyConstraint();
    }
  return(*this);
}

ION_FILTERED_ITERATOR_T& ION_FILTERED_ITERATOR_T::operator++(int c) {
  for (int i=0;i<c;i++) {
    if (current != ion_series -> end()) {
      ++current;
      satisfyConstraint();
    }
  }
  return(*this);
}

ION_T*& ION_FILTERED_ITERATOR_T::operator*() {
  return *current;
}

ION_T* ION_FILTERED_ITERATOR_T::operator->() {
  return *current;
}

/*
 * if the current ion doesn't satisfy the constraint
 * find the next one until we reached the end of the 
 * list.  If we reach the end of the list, set current
 * to the ion_series's end so that we can use end
 * in all of the iterating loops.
 */
void ION_FILTERED_ITERATOR_T::satisfyConstraint() {
  if (current == ion_series -> end()) 
    return;

  while (current != end_iter && !constraint -> is_satisfied(*current)) {
    ++current;    
  }

  if (current == end_iter || current == ion_series -> end()) {
    current = ion_series -> end();
  }
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
