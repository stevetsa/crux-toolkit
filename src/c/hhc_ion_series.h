#ifndef HHC_ION_SERIES_H
#define HHC_ION_SERIES_H
#include "hhc.h"
#include "scorer.h"
#include "ion_series.h"
#include "spectrum.h"
#include "spectrum_collection.h"

#define bin_width_mono 1.0005079

using namespace std;

class LinkedIonSeries {
  public:
    LinkedIonSeries() {}
    LinkedIonSeries(char* links, int charge, FLOAT_T linker_mass);
    LinkedIonSeries(char* sequenceA, char* sequenceB, int posA, int posB, int charge, FLOAT_T linker_mass);
    LinkedIonSeries(char* sequenceA, char* sequenceB, char* links, int charge, FLOAT_T linker_mass);
    int charge() { return charge_; }
    void set_charge(int charge) { charge_ = charge; }
    std::vector<LinkedPeptide>& ions() { return all_ions; }
    void add_linked_ions(LinkedPeptide& linked_peptide);
    //void add_linked_ions_double(char* sequenceA, char* sequenceB);
    //void add_linked_ions_self(char* sequence);
    //void add_linked_ions_self_loop(char* sequence);
    //void add_linked_ions_self_dead_end(char* sequence);
    int size() { return all_ions.size(); }
    void print();
    //int num_peptides;
  private:
    FLOAT_T linker_mass_;
    std::map<char, std::set<char> > bond_map;
    int charge_; 
    // a list of all the ions 
    std::vector<LinkedPeptide> all_ions;

};

LinkedIonSeries::LinkedIonSeries(char* links, int charge, FLOAT_T linker_mass) {
  charge_ = charge;
  linker_mass_ - linker_mass;
  //num_peptides = 0;
  std::string bonds_string = std::string(links);
  for (int i = 0; i < bonds_string.length() - 2; i += 4) {
     bond_map[bonds_string[i]].insert(bonds_string[i+2]);
     bond_map[bonds_string[i+2]].insert(bonds_string[i]);
  }
}



void LinkedIonSeries::print() {
  sort(all_ions.begin(), all_ions.end());
  for (vector<LinkedPeptide>::iterator it = all_ions.begin(); it != all_ions.end(); ++it) {
    cout << it->get_mz() << " " << *it << endl;
  }
}
void LinkedIonSeries::add_linked_ions(LinkedPeptide& linked_peptide) {
  std::vector<std::pair<LinkedPeptide, LinkedPeptide> > fragments;
  //++num_peptides;
  linked_peptide.split(fragments);
  for (vector<pair<LinkedPeptide, LinkedPeptide> >::iterator it = fragments.begin(); it != fragments.end(); ++it) {
    if (it->first.charge() != 0) {
      it->first.set_type(B_ION); 
      it->first.calculate_mass();
      it->first.get_mz();
      //cout << it->first.get_mz() << " B " << it->first << endl;
      all_ions.push_back(it->first);
    }
    if (it->second.charge() != 0) {
      it->second.set_type(Y_ION); 
      it->second.calculate_mass();
      it->second.get_mz();
      //cout << it->second.get_mz() << " Y " << it->second << endl;
      all_ions.push_back(it->second);
    }
  }
}

LinkedIonSeries::LinkedIonSeries(char* sequenceA, char* sequenceB, int posA, int posB, int charge, FLOAT_T linker_mass) {
  charge_ = charge;
  linker_mass_ = linker_mass;
  Peptide p = Peptide(sequenceA);
  //LinkedPeptide lp = LinkedPeptide(p, charge, linker_mass);
  LinkedPeptide lp = LinkedPeptide(sequenceA, sequenceB, posA, posB, linker_mass, charge);
  add_linked_ions(lp);
}


LinkedIonSeries::LinkedIonSeries(char* sequenceA, char* sequenceB, char* links, int charge, FLOAT_T linker_mass) {
    charge_ = charge;
    linker_mass_ = linker_mass;
    std::string bonds_string = std::string(links);
    for (int i = 0; i < bonds_string.length() - 2; i += 4) {
       bond_map[bonds_string[i]].insert(bonds_string[i+2]);
    }

  std::string alpha_sequence = string(sequenceA);
  std::string beta_sequence = string(sequenceB);
  for (int i = 0; i < alpha_sequence.length(); ++i) {
     std::map<char, std::set<char> >::iterator char_it = bond_map.find(alpha_sequence[i]);
     if (char_it != bond_map.end()) {
       for (int j = 0; j < beta_sequence.length(); ++j) {
         if (char_it->second.find(beta_sequence[j]) != char_it->second.end()) {
           //cout << alpha_sequence << " " << i << " " << beta_sequence << " " << j << endl;
	   LinkedPeptide lp = LinkedPeptide(sequenceA,sequenceB,i,j,linker_mass,charge);
	   //cout << "new linked peptide: " << lp << endl;
           add_linked_ions(lp);
          }
         }
        }
      }
   for (int i = 0; i < beta_sequence.length(); ++i) {
     std::map<char, set<char> >::iterator char_it = bond_map.find(beta_sequence[i]);
     if (char_it != bond_map.end()) {
       for (int j = 0; j < alpha_sequence.length(); ++j) {
         if (char_it->second.find(alpha_sequence[j]) != char_it->second.end()) {
	   LinkedPeptide lp = LinkedPeptide(sequenceB,sequenceA,i,j,linker_mass,charge);
           add_linked_ions(lp);
         }
       }	
     }
    }
    
}


//void print_observed_spectrum(SPECTRUM_T* spectrum);

void print_spectrums(FLOAT_T* theoretical, SPECTRUM_T* spectrum, FLOAT_T min_mz_float, FLOAT_T max_mz_float, int scale);

bool hhc_create_intensity_array_theoretical(
  SCORER_T* scorer,        ///< the scorer object -in/out
  LinkedIonSeries& ion_series,
  FLOAT_T* theoretical       ///< the empty theoretical spectrum -out
  );

FLOAT_T hhc_score_spectrum_v_ion_series(
  SCORER_T* scorer,        ///< the scorer object -in
  SPECTRUM_T* spectrum,    ///< the spectrum to score -in
  LinkedIonSeries& ion_series ///< the ion series to score against the spectrum -in
  );

FLOAT_T hhc_gen_score_xcorr(
  SCORER_T* scorer,        ///< the scorer object -in
  SPECTRUM_T* spectrum,    ///< the spectrum to score -in
  LinkedIonSeries& ion_series ///< the ion series to score against the spectrum -in
  );

FLOAT_T hhc_gen_score_xcorr(
  SCORER_T* scorer,        ///< the scorer object -in
  SPECTRUM_T* spectrum,    ///< the spectrum to score -in
  LinkedIonSeries& ion_series ///< the ion series to score against the spectrum -in
  )
{
  FLOAT_T final_score = 0;
  FLOAT_T* theoretical = NULL;

  // initialize the scorer before scoring if necessary
  // preprocess the observed spectrum in scorer
  if(!scorer->initialized){
    // create intensity array for observed spectrum, if already not been done
    //if(!create_intensity_array_xcorr(spectrum, scorer, get_ion_series_charge(ion_series))){
    
    if (!create_intensity_array_xcorr(spectrum, scorer, ion_series.charge())) {
      carp(CARP_FATAL, "failed to produce XCORR");
    }
  }
  
  // create theoretical array
  theoretical = (FLOAT_T*)mycalloc(scorer->sp_max_mz, sizeof(FLOAT_T));
  
  // create intensity array for theoretical spectrum 
  //if(!create_intensity_array_theoretical(scorer, ion_series, theoretical)){
  if (!hhc_create_intensity_array_theoretical(scorer, ion_series, theoretical)) {
    carp(CARP_ERROR, "failed to create theoretical spectrum for Xcorr");
    return FALSE;
  }
  
  // do cross correlation between observed spectrum(in scorer) and theoretical spectrum.
  // use the two intensity arrays that were created
  final_score = cross_correlation(scorer, theoretical);
  //print_spectrums(theoretical, spectrum, spectrum->max_peak_mz, 1);
  //print_spectrums(theoretical, spectrum, 300, 750, 1);

  // free theoretical spectrum
  free(theoretical);

  // return score
  return final_score;
}

FLOAT_T hhc_score_spectrum_v_ion_series(
  SCORER_T* scorer,        ///< the scorer object -in
  SPECTRUM_T* spectrum,    ///< the spectrum to score -in
  LinkedIonSeries& ion_series ///< the ion series to score against the spectrum -in
  )
{
  FLOAT_T final_score = 0;
  // if score type equals SP
  if(scorer->type == SP){
    //final_score = gen_score_sp(scorer, spectrum, ion_series);
  }
  else if(scorer->type == XCORR){
    final_score = hhc_gen_score_xcorr(scorer, spectrum, ion_series);
  }
  // FIXME, later add different score types...
  else{
    carp(CARP_ERROR, "no scoring method availiable for the scorers' score type");
  }
  return final_score;
}

bool hhc_create_intensity_array_theoretical(
  SCORER_T* scorer,        ///< the scorer object -in/out
  LinkedIonSeries& ion_series,
  FLOAT_T* theoretical       ///< the empty theoretical spectrum -out
  )
{
  //ION_T* ion = NULL;
  int ion_charge = 0;
  ION_TYPE_T ion_type;
  int intensity_array_idx = 0;
  FLOAT_T bin_width = bin_width_mono;
  vector<LinkedPeptide>& ions = ion_series.ions();
  // while there are ion's in ion iterator, add matched observed peak intensity
  for (vector<LinkedPeptide>::iterator ion = ions.begin(); ion != ions.end(); ++ion) {
  //while(ion_iterator_has_next(ion_iterator)){
    //ion = ion_iterator_next(ion_iterator);
    //ion->calculate_mass();
    intensity_array_idx = (int)(ion->get_mz() / bin_width + 0.5);
    //cout << "index " << intensity_array_idx << endl;
    ion_type = ion->type();
    ion_charge = ion->charge();
    //cout << "m/z: " << ion->get_mz() << " charge: " << ion->charge() << endl;
    // skip ions that are located beyond max mz limit
    if(intensity_array_idx >= scorer->sp_max_mz){
      continue;
    }
  //if (ion->type() == B_ION) cout << "b-ion" << endl;
  //if (ion->type() == Y_ION) cout << "y-ion" << endl;
  //cout << *ion << endl;
  //if (ion->type() == B_ION) 
    //cout << "bion" << endl; else cout << "yion" << endl;
/*
    if  (ion_type == B_ION){
      if(ion_is_modified(ion)){
        carp(CARP_INFO, "idx: %d, adding ion type: MOD-%s",  intensity_array_idx,  "B");
      }
      else{
        carp(CARP_INFO, "idx: %d, adding ion type: %s",  intensity_array_idx,  "B");
      }
    }
    else if(ion_type == Y_ION){
      if(ion_is_modified(ion)){
        carp(CARP_INFO, "idx: %d, adding ion type: MOD-%s",  intensity_array_idx,  "Y");
      }
      else{
        carp(CARP_INFO, "idx: %d, adding ion type: %s",  intensity_array_idx,  "Y");
      }
    }
    else{
      carp(CARP_INFO, "idx: %d, adding ion type: %s",  intensity_array_idx,  "A");
    }
*/

    // is it B, Y ion?

    // neutral loss peak?
    // Add peaks of intensity 50.0 for B, Y type ions. 
    // In addition, add peaks of intensity of 25.0 to +/- 1 m/z flanking each B, Y ion.
    // Skip ions that are located beyond max mz limit
    if((intensity_array_idx)< scorer->sp_max_mz){
      //if (ion->type() == Y_ION)
      //add_intensity(theoretical, intensity_array_idx, 51);
      //else 
      add_intensity(theoretical, intensity_array_idx, 50);
      add_intensity(theoretical, intensity_array_idx - 1, 25);
    }
    if((intensity_array_idx + 1)< scorer->sp_max_mz){
      add_intensity(theoretical, intensity_array_idx + 1, 25);
    }

    // add neutral loss of water and NH3
    // mass_z + (modification_masses[(int)ion_modification]/(FLOAT_T)charge) * modification_count;  


    if(ion_type == B_ION){
      int h2o_array_idx = (int)((ion->get_mz() - (MASS_H2O_MONO/ion->charge()) ) / bin_width + 0.5);
      add_intensity(theoretical, h2o_array_idx, 10);
    }

    int nh3_array_idx = (int)((ion->get_mz() -  (MASS_NH3_MONO/ion->charge())) / bin_width + 0.5);
    add_intensity(theoretical, nh3_array_idx, 10);        
  }
  return true;
}


#endif
