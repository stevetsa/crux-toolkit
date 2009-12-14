#include "xhhc.h"
#include "xhhc_ion_series.h"
#include "xhhc_scorer.h"

extern "C" {
#include "objects.h"
#include "scorer.h"
#include "spectrum_collection.h"
}


#include <math.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#define bin_width_mono 1.0005079


#define NUM_ARGUMENTS 8
#define NUM_OPTIONS 4


double get_concat_score(char* peptideA, char* peptideB, int link_site, int charge, SPECTRUM_T* spectrum);
void print_spectrum(SPECTRUM_T* spectrum, LinkedIonSeries& ion_series);
int main(int argc, char** argv){

  /* Verbosity level for set-up/command line reading */
  set_verbosity_level(CARP_ERROR);
  
  /* Define optional command line arguments */
  int num_options = NUM_OPTIONS;
  const char* option_list[NUM_OPTIONS] = {
    "verbosity",
    "version",
    "xcorr-use-flanks",
    "xlink-score-method"
  };

  

  /* Define required command line arguments */
  int num_arguments = NUM_ARGUMENTS;
  const char* argument_list[NUM_ARGUMENTS] = {"peptide A",
                                              "peptide B",
					      "pos A",
					      "pos B",
					      "link mass",
					      "charge state",
					      "scan number",
					      "ms2 file"};

  
  /* for debugging of parameter processing */
  //set_verbosity_level( CARP_DETAILED_DEBUG );
  set_verbosity_level( CARP_ERROR );
  
  /* Set default values for parameters in parameter.c */
  initialize_parameters();

  /* Define optional and required command line arguments */
  select_cmd_line_options( option_list, num_options );
  select_cmd_line_arguments( argument_list, num_arguments);

  /* Parse the command line, including the optional params file */
  /* does sytnax, type, bounds checking and dies if neccessessary */
  parse_cmd_line_into_params_hash(argc, argv, "xlink-score-peptide-spectrum");

  /* Set verbosity */
  set_verbosity_level(get_int_parameter("verbosity"));

  /* Get Arguments */
  char* peptideA = get_string_parameter("peptide A");
  char* peptideB = get_string_parameter("peptide B");
  
  int posA     = get_int_parameter("pos A");
  int posB     = get_int_parameter("pos B");
  int charge   = get_int_parameter("charge state"); 
  int scan_num = get_int_parameter("scan number"); 

  char* ms2_file = get_string_parameter("ms2 file");

  LinkedPeptide::linker_mass = get_double_parameter("link mass");
 
  // create new ion series
  
  // a single peptide linked to itself
  if (strcmp(peptideB, "NULL") == 0) {
    cout << "B is null" << endl; 
    peptideB = NULL;
  }

  // read ms2 file
  SPECTRUM_COLLECTION_T* collection = new_spectrum_collection(ms2_file);
  SPECTRUM_T* spectrum = allocate_spectrum();
  //cout << "lp " << lp << endl; 
  // search for spectrum with correct scan number
  if(!get_spectrum_collection_spectrum(collection, scan_num, spectrum)){
    carp(CARP_ERROR, "failed to find spectrum with  scan_num: %d", scan_num);
    free_spectrum_collection(collection);
    free_spectrum(spectrum);
    exit(1);
  }

  //created linked peptide.
  LinkedPeptide lp = LinkedPeptide(peptideA, peptideB, posA, posB, charge);
  Scorer xhhc_scorer;
  xhhc_scorer.set_print(false);

  string scoremethod(get_string_parameter("xlink-score-method"));

  if (scoremethod=="composite") {

    LinkedIonSeries ion_series;

    //cout << lp << endl;
    
    ion_series.add_linked_ions(lp);
     
    double score = xhhc_scorer.score_spectrum_vs_series(spectrum, ion_series);

    cout <<score<<endl;

    bool do_print_spectra = true;
    if (do_print_spectra) {
      print_spectrum(spectrum, ion_series);
    }
  } else if (scoremethod=="modification") {
    
    LinkedIonSeries ion_seriesA;
    ion_seriesA.add_linked_ions(lp, 1);
    double scoreA = xhhc_scorer.score_spectrum_vs_series(spectrum, ion_seriesA);
    
    LinkedIonSeries ion_seriesB;
    ion_seriesB.add_linked_ions(lp, 2);

    

    double scoreB = xhhc_scorer.score_spectrum_vs_series(spectrum, ion_seriesB);

    if (scoreA > scoreB)
      cout << scoreA << "\t" << scoreB << endl;
    else
      cout << scoreB << "\t" << scoreA << endl;

  } else if (scoremethod=="concatenation") {


    vector<double> scores;
    double score1 = get_concat_score(peptideA, peptideB, posA, charge, spectrum);
    scores.push_back(score1);

    double score2 = get_concat_score(peptideB, peptideA, posB, charge, spectrum);
    scores.push_back(score2);


    int lengthA = string(peptideA).length();
    int lengthB = string(peptideB).length();

    double score3 = get_concat_score(peptideA, peptideB, lengthA + posB, charge, spectrum);
    scores.push_back(score3);

    double score4 = get_concat_score(peptideB, peptideA, lengthB + posA, charge, spectrum);
    scores.push_back(score4);

    sort(scores.begin(), scores.end(), less<double>());
    cout <<scores[0];
    for (int i=1;i<4;i++)
      {
	cout <<"\t"<<scores[i];
      }

    cout << endl;
  }
  else {
    carp(CARP_ERROR,"Unknown method");
  }
  // free heap
  free_spectrum_collection(collection);
  free_spectrum(spectrum);
}


double get_concat_score(char* peptideA, char* peptideB, int link_site, int charge, SPECTRUM_T* spectrum) {
  string lpeptide = string(peptideA) + string(peptideB); 
  
  ION_CONSTRAINT_T* ion_constraint = new_ion_constraint_smart(XCORR, charge);
  
  ION_SERIES_T* ion_series = new_ion_series(lpeptide.c_str(), charge, ion_constraint);
  
  predict_ions(ion_series);
  
  //modify ions.

  ION_ITERATOR_T* ion_iterator = new_ion_iterator(ion_series);
  
  
  //int pepA_begin = 0;
  int pepB_begin = string(peptideA).length();
  int llength = lpeptide.length();
  
  while(ion_iterator_has_next(ion_iterator)){
    ION_T* ion = ion_iterator_next(ion_iterator);
      //check to see if if is the cterm of 1st peptide.
      int ion_charge = get_ion_charge(ion);
      int cleavage_idx = get_ion_cleavage_idx(ion);
      ION_TYPE_T ion_type = get_ion_type(ion);

      //if contains cterm of 1st peptide, modify by -OH 
      
      carp(CARP_DEBUG,"====================");
      if (ion_type == B_ION) {
	carp(CARP_DEBUG,"B-ion");
	carp(CARP_DEBUG,"%s",lpeptide.substr(0,cleavage_idx).c_str());
      } else if (ion_type == Y_ION) {
	carp(CARP_DEBUG,"Y-ion");
	carp(CARP_DEBUG,"%s",lpeptide.substr(llength-cleavage_idx,llength).c_str());
      }
      else continue;

      carp(CARP_DEBUG,"cleavage idx:%d",cleavage_idx);
      //print_ion(ion, stdout);

      bool cterm_1st = false;
      if (ion_type == B_ION) {
	carp(CARP_DEBUG,"B-Ion");
	if (cleavage_idx >= pepB_begin) {
	  cterm_1st = true;
	}
      } else if (ion_type == Y_ION) {
	carp(CARP_DEBUG,"Y-ion");
	if (cleavage_idx > (llength - pepB_begin)) {
	  cterm_1st = true;
	}
      }

      bool nterm_2nd = false;
      if (ion_type == B_ION) {
	if (cleavage_idx > pepB_begin) {
	  nterm_2nd = true;
	}
      } else if (ion_type == Y_ION) {
	if (cleavage_idx >= (llength-pepB_begin)) {
	  nterm_2nd = true;
	}
      }

      bool has_link_site = false;
     if (ion_type == B_ION) {
	if (cleavage_idx > link_site) {
	  has_link_site = true;
	}
      } else if (ion_type == Y_ION) {
       if (cleavage_idx >= (llength- link_site)) {
	  has_link_site = true;
	}
      }
      
     carp(CARP_DEBUG,"cterm:%d",cterm_1st);
     carp(CARP_DEBUG,"nterm:%d",nterm_2nd);
     carp(CARP_DEBUG,"has site:%d",has_link_site);
      

      //if it contains the cterm of the 1st peptide, modify by -OH
      if (cterm_1st) {
	FLOAT_T old_mass = (get_ion_mass_z(ion) - MASS_H_MONO) * (FLOAT_T)ion_charge;
	FLOAT_T new_mass = old_mass + MASS_H2O_MONO - MASS_H_MONO;
	FLOAT_T new_mz = (new_mass + (FLOAT_T)ion_charge) / (FLOAT_T)ion_charge;
	set_ion_mass_z(ion, new_mz);
      }
      //if contains the nterm of 2nd peptide, modify by -H
      if (nterm_2nd) {
	FLOAT_T old_mass = (get_ion_mass_z(ion) - MASS_H_MONO) * (FLOAT_T)ion_charge;
	FLOAT_T new_mass = old_mass + MASS_H_MONO;
	FLOAT_T new_mz = (new_mass + (FLOAT_T)ion_charge) / (FLOAT_T)ion_charge;
	set_ion_mass_z(ion, new_mz);
      }
      //if contains the link site, modify by link mass.
      if (has_link_site) {
	FLOAT_T old_mass = (get_ion_mass_z(ion) - MASS_H_MONO) * (FLOAT_T)ion_charge;
	FLOAT_T new_mass = old_mass + LinkedPeptide::linker_mass;
	FLOAT_T new_mz = (new_mass + (FLOAT_T)ion_charge) / (FLOAT_T)ion_charge;
	set_ion_mass_z(ion, new_mz);
      }
    

    

      //print_ion(ion, stdout);
    
    }

    SCORER_T* scorer = new_scorer(XCORR); 

    // calculate the score
    FLOAT_T score = score_spectrum_v_ion_series(scorer, spectrum, ion_series);
    return score;



}

void print_spectrum(SPECTRUM_T* spectrum, LinkedIonSeries& ion_series) {

      int total_by_ions = ion_series.get_total_by_ions();
      int matched_by_ions = Scorer::get_matched_by_ions(spectrum, ion_series);
      FLOAT_T frac_by_ions = (double)matched_by_ions / (double) total_by_ions;

      carp(CARP_DEBUG,"total:%d",total_by_ions);
      carp(CARP_DEBUG,"matched:%d",matched_by_ions);
      carp(CARP_DEBUG,"frac:%d",frac_by_ions);
      carp(CARP_DEBUG,"npeaks:",get_spectrum_num_peaks(spectrum));
      FLOAT_T bin_width = bin_width_mono;
      vector<LinkedPeptide>& ions = ion_series.ions();
      
      map<PEAK_T*, bool> matched;

      for (vector<LinkedPeptide>::iterator ion = ions.begin();
	   ion != ions.end(); 
	   ++ion) {
	if (ion -> get_mz(MONO) >= 400 && ion -> get_mz(MONO) <= 1200) {
	  if (ion -> type() == B_ION || ion -> type() == Y_ION) {
	    PEAK_T* peak = get_nearest_peak(spectrum, ion -> get_mz(MONO), bin_width);
	    if (peak != NULL) {
	      matched[peak] = true;
	    }
	  }
	}
      }

      //now print the spectrum and whether or not it has been matched.
      
      ofstream fout("ion_match.out");

      PEAK_ITERATOR_T* peak_iter = new_peak_iterator(spectrum);
      while (peak_iterator_has_next(peak_iter)) {
	PEAK_T* peak = peak_iterator_next(peak_iter);

	if (get_peak_location(peak) >= 400 && get_peak_location(peak) <= 1200) {

	  fout << get_peak_location(peak) << "\t";
	  fout << get_peak_intensity(peak) << "\t";
	  if (matched.find(peak) == matched.end())
	    fout << "0";
	  else
	    fout << "1";
	  fout<<endl;
	}
      }

      fout.close();
      free_peak_iterator(peak_iter);



}
