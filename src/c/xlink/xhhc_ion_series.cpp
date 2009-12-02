#include "xhhc_ion_series.h"
#include "xhhc_scorer.h"

#include <iostream>

// main constructor
LinkedIonSeries::LinkedIonSeries(char* links, int charge) {
  charge_ = charge;
  string bonds_string = string(links);
  int stop = (int) bonds_string.length() - 2;
  for (int i = 0; i < stop; i += 4) {
     bond_map[bonds_string[i]].insert(bonds_string[i+2]);
     bond_map[bonds_string[i+2]].insert(bonds_string[i]);
  }
}

// prints out tab delimited information about the ion series
void LinkedIonSeries::print() {
  sort(all_ions.begin(), all_ions.end());
  cout << "m/z\ttype\tion" << endl;
  string ion_type;
  for (vector<LinkedPeptide>::iterator ion = all_ions.begin(); ion != all_ions.end(); ++ion) {
    if (ion->type() == B_ION) 
      ion_type = "B_ION";
    else 
      ion_type = "Y_ION";
    cout << ion->get_mz(get_mass_type_parameter("fragment-mass")) << "\t" << ion_type << "\t" << *ion << endl;
  }
}

// cleaves linked_peptide at all positions, adding b and y ions
void LinkedIonSeries::add_linked_ions(LinkedPeptide& linked_peptide) {
  if (charge_ == 0) charge_ = linked_peptide.charge();
  vector<pair<LinkedPeptide, LinkedPeptide> > fragments;
  // split the precursor at every cleavage site
  linked_peptide.split(fragments);
  for (vector<pair<LinkedPeptide, LinkedPeptide> >::iterator ion_pair = fragments.begin(); ion_pair != fragments.end(); ++ion_pair) {
    // if b-ion and not a neutral loss
    if (ion_pair->first.charge() != 0) {
      ion_pair->first.set_type(B_ION); 
      //ion_pair->first.calculate_mass();
      //cout << ion_pair->first.get_mz() << " B " << ion_pair->first << endl;
      all_ions.push_back(ion_pair->first);
    }
    // if y-ion and not a neutral loss
    if (ion_pair->second.charge() != 0) {
      ion_pair->second.set_type(Y_ION); 
      //ion_pair->second.calculate_mass();
      //cout << ion_pair->second.get_mz() << " Y " << ion_pair->second << endl;
      all_ions.push_back(ion_pair->second);
    }
  }
}



/***************************************
 *CRUX OVERRIDES
 ***************************************/


/**
 * \brief Creates an array in which element i is the sum of the masses
 * of amino acids 0 to (i-1).  At i=0 is stored the length of the
 * peptide.  
 * \returns an array of ion masses for all sub sequences
 */
FLOAT_T* hhc_create_ion_mass_matrix(
  //char* peptide, ///< The peptide for this ion series. -in
  MODIFIED_AA_T* modified_seq, ///< the sequence
  MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
  int peptide_length, ///< the length of the peptide
  FLOAT_T linker_mass,
  int linker_site
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
  // for open modification cross linking
  if (linker_site == 0) mass_matrix[ion_idx] += linker_mass; 
  //++ion_idx;
  //for(; ion_idx <= peptide_length; ++ion_idx){
  for(ion_idx = 2; ion_idx <= peptide_length; ++ion_idx){
    mass_matrix[ion_idx] = mass_matrix[ion_idx-1] + 
      get_mass_mod_amino_acid(modified_seq[ion_idx-1], mass_type);
    if (linker_site == ion_idx-1)
      mass_matrix[ion_idx] += linker_mass; 
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
 * \brief The engine of ion series. Predicts all the ions from the
 * peptide that meet the ion constraint. All predicted ions are stored
 * in the ion_series as ion objects. 
 */
void hhc_predict_ions(
  ION_SERIES_T* ion_series, ///< the ion series to predict ions for -in
  FLOAT_T linker_mass,
  int linker_site
  )
{


  if(get_ion_series_is_predicted(ion_series)){
    carp(CARP_WARNING, "The ion series has already been predicted and added");
    return;
  }


  ION_CONSTRAINT_T* constraint = get_ion_series_ion_constraint(ion_series);
  
  // create a mass matrix
  FLOAT_T* mass_matrix = 
    hhc_create_ion_mass_matrix(get_ion_series_modified_aa_seq(ion_series), 
			       get_ion_constraint_mass_type(constraint), 
			       get_ion_series_peptide_length(ion_series), 
			       linker_mass, linker_site);  
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
  scan_for_aa_for_neutral_loss(ion_series);
  
  // generate ions without any modifications
  if(!generate_ions_no_modification(ion_series, mass_matrix)){
    carp(CARP_FATAL, "failed to generate ions, no modifications");
  }

  // create modification ions?
  if(get_ion_constraint_use_neutral_losses(constraint)){
    
    // generate ions with nh3 modification
    if(abs(get_ion_constraint_modification(constraint, NH3)) > 0){
      if(!generate_ions(ion_series, NH3)){
        carp(CARP_FATAL, "failed to generate ions, NH3 modifications");
      }
    }
    
    // generate ions with h2o modification
    if(abs(get_ion_constraint_modification(constraint, H2O)) > 0){
      if(!generate_ions(ion_series, H2O)){
        carp(CARP_FATAL, "failed to generate ions, H2O modifications");
      }
    }

    // generate ions with isotope modification
    if(get_ion_constraint_modification(constraint,ISOTOPE) > 0){
      if(!generate_ions(ion_series, ISOTOPE)){
        carp(CARP_FATAL, "failed to generate ions, ISOTOPE modifications");
      }
    }

    // generate ions with flank modification
    if(get_ion_constraint_modification(constraint,FLANK) > 0){
      if(!generate_ions_flank(ion_series)){
        carp(CARP_FATAL, "failed to generate ions, FLANK modifications");
      }
    }
    
    // add more modifications here

  }
  
  // ion series now been predicted
  set_ion_series_is_predicted(ion_series, TRUE);

  // free mass matrix
  free(mass_matrix);
}
