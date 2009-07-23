#include <math.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include "hhc_ion_series.h"
#include "objects.h"

#define bin_width_mono 1.0005079


int noflanks = 0;
int normalize = 0;

int main(int argc, char** argv){

  // required variables
  char* ms2_file = NULL;
  int scan_num = 0;
  char* positions = NULL;
  char* sequenceA = NULL;
  char* sequenceB = NULL;
  char* links = NULL;
  // optional variables
  char* linker_mass_string = "0";
  int charge = 2;
  //char* type = "xcorr";
  char* parameter_file = "crux.params";
  int  verbosity = CARP_ERROR;
  // parsing variables
  int result = 0;
  char * error_message; 
  /* Define optional command line arguments */ 

  parse_arguments_set_opt(
    "charge", 
    "The peptide charge. 1|2|3",
    (void *) &charge, 
    INT_ARG);
 /* 
  parse_arguments_set_opt(
    "score-type", 
    "The type of scoring function to use. sp | xcorr",
    (void *) &type, 
    STRING_ARG);
*/
  parse_arguments_set_opt(
    "parameter-file",
    "The crux parameter file to parse parameter from.",
    (void *) &parameter_file,
    STRING_ARG);

  parse_arguments_set_opt(
    "link-at", 
    "Specific link positions on peptide A and B ex. 4,0",
    (void *) &positions, 
    STRING_ARG);

  parse_arguments_set_opt(
    "links", 
    "comma delimited pair of amino acid link sites, ex. A:K,A:D", 
    (void *) &links, 
    STRING_ARG);

  parse_arguments_set_opt(
    "linker-mass", 
    "linker mass and linker modifications. Default 0.",
    (void *) &linker_mass_string, 
    STRING_ARG);

  parse_arguments_set_opt(
    "no-flanks", 
    "graph theoretical spectrum without flanks",
    (void *) &noflanks, 
    INT_ARG);
  
  parse_arguments_set_opt(
    "normalize", 
    "normalize theoretical spectrum intensities",
    (void *) &normalize, 
    INT_ARG);

  parse_arguments_set_opt(
    "verbosity", 
    "Specify the verbosity of the current processes from 0-100.",
    (void *) &verbosity, 
    INT_ARG);
  
 /* Define required command line arguments */
  parse_arguments_set_req(
    "peptide-sequence-alpha", 
    "The first peptide sequence.", 
    (void *) &sequenceA, 
    STRING_ARG);

  parse_arguments_set_req(
    "peptide-sequence-beta", 
    "The second peptide sequence.", 
    (void *) &sequenceB, 
    STRING_ARG);
  
  parse_arguments_set_req(
    "scan-number", 
    "The scan number for the MS-MS spectrum to extract from the ms2 file. This is an integer in the range [1, 100000], and uniquely identifies a particular MS-MS spectrum within an .ms2 file.",
    (void *) &scan_num, INT_ARG);

  parse_arguments_set_req(
    "ms2-filename", 
    "A file containing multiple MS-MS spectra in .ms2 format.",
    (void *) &ms2_file,
    STRING_ARG);

  /* Parse the command line */
  if (parse_arguments(argc, argv, 0)) {
    // parsed arguments
    int peptide_charge = 1;
    SCORER_TYPE_T score_type = XCORR; 
    
    SPECTRUM_T* spectrum = NULL;
    SPECTRUM_COLLECTION_T* collection = NULL;
    SCORER_T* scorer = NULL;
    FLOAT_T score = 0;
   initialize_parameters();    
    FLOAT_T wrong = calc_sequence_mass(sequenceB, MONO);
    FLOAT_T right = calc_mod_sequence_mass(sequenceB, MONO);    

   cout << "mass: " << wrong << " mono mass: " << right << endl;
    if (positions == NULL && links == NULL) {
      carp(CARP_FATAL, "must specify either link amino acids or specific link position");
    }
    // set verbosity
    if(CARP_FATAL <= verbosity && verbosity <= CARP_MAX){
      set_verbosity_level(verbosity);
    }
    else{
      carp(CARP_FATAL, "verbosity level must be between 0-100");
    }

    // set verbosity
    //set_verbosity_level(verbosity);
    
    //peptide_charge = get_int_parameter("charge");
    
    if( charge < 1 || charge > 5){
      carp(CARP_FATAL, "peptide charge must be between 1 and 5.");
    }

    // check peptide sequence
    if(!valid_peptide_sequence(sequenceA)) {
      carp(CARP_ERROR, "%s not a valid sequence", sequenceA);
    }

    if(!valid_peptide_sequence(sequenceB) && strcmp(sequenceB, "NULL") != 0) {
      carp(CARP_ERROR, "%s not a valid sequence", sequenceB);
    }
    
    // score type
    /*
    if(strcmp(get_string_parameter_pointer("score-type"), "sp")== 0){
      score_type = SP;
    }
    else if(strcmp(get_string_parameter_pointer("score-type"), "xcorr")== 0){
      score_type = XCORR;
    }
    else{
      carp(CARP_ERROR, "The type of scoring function to use. sp | xcorr");
      //wrong_command(type, "The type of scoring function to use. sp | xcorr");
    }
     */ 
    // parameters are now confirmed, can't be changed
    //parameters_confirmed();
    
    // set ion constraint to sequest settings
    //ION_CONSTRAINT_T* ion_constraint = NULL;
  /*  
    if(score_type == SP){
      //ion_constraint = new_ion_constraint_sequest_sp(peptide_charge);  
      // create new scorer
      scorer = new_scorer(SP);  
    }
    else if(score_type == XCORR){
      //ion_constraint = new_ion_constraint_sequest_xcorr(peptide_charge);  
      scorer = new_scorer(XCORR);  
    }
  */ scorer = new_scorer(XCORR);
    // create new ion series
    LinkedIonSeries ion_series; 
    FLOAT_T linker_mass = atof(linker_mass_string);
    // a single peptide linked to itself
    if (strcmp(sequenceB, "NULL") == 0) {
      cout << "B is null" << endl; 
      sequenceB = NULL;
    }
    if (positions != NULL) {
      char * comma = strchr(positions, ',');
      *comma++ = '\0';      
      int posA =  atoi(positions);
      int posB = atoi(comma);
      ion_series = LinkedIonSeries(sequenceA, sequenceB, posA, posB, charge, linker_mass);
    } else {
      ion_series = LinkedIonSeries(sequenceA, sequenceB, links, charge, linker_mass);
    }
    //cout << ion_series.size();
    carp(CARP_DEBUG, "number of ions: %d\n", ion_series.size());
    //ion_series = new_ion_series(peptide_sequence, peptide_charge, ion_constraint);
    ion_series.print(); 
   // now predict ions
    //predict_ions(ion_series);
      
   // read ms2 file
    collection = new_spectrum_collection(ms2_file);
    spectrum = allocate_spectrum();
    
    // search for spectrum with correct scan number
    if(!get_spectrum_collection_spectrum(collection, scan_num, spectrum)){
      carp(CARP_ERROR, "failed to find spectrum with  scan_num: %d", scan_num);
      //free_ion_constraint(ion_constraint);
      //free_ion_series(ion_series);
      free_spectrum_collection(collection);
      free_spectrum(spectrum);
      exit(1);
    }
    
    // calculates the Sp score
    //score = score_spectrum_v_ion_series(scorer, spectrum, ion_series);
    score = hhc_score_spectrum_v_ion_series(scorer, spectrum, ion_series);
    // print the Sp score
    if(score_type == SP){
      printf("Sp score is: %.2f\n", score);
    }
    else if(score_type == XCORR){
      printf("Xcorr score is: %.2f\n", score);
    }
    else{
      carp(CARP_ERROR, "invalid score type for the application");
    }
      
   // free heap
   free_scorer(scorer);
   //free_ion_constraint(ion_constraint);
   //free_ion_series(ion_series);
   free_spectrum_collection(collection);
   free_spectrum(spectrum);
   free_parameters();
 }
 else{
   char* usage = parse_arguments_get_usage("score_peptide_spectrum");
   result = parse_arguments_get_error(&error_message);
   fprintf(stderr, "Error in command line. Error # %d\n", result);
   fprintf(stderr, "%s\n", error_message);
   fprintf(stderr, "%s", usage);
   free(usage);
 }
 exit(0);
}

//creates three files for spectacle.pl: spectrums.out, theoretical.out, observed.out
void print_spectrums(FLOAT_T* theoretical, SPECTRUM_T* spectrum, FLOAT_T min_mz_float, FLOAT_T max_mz_float, int scale) {
   
  ofstream theoretical_file;
  ofstream observed_file;
  ofstream spectrums_file;
  theoretical_file.open("theoretical.out");
  observed_file.open("observed.out");
  spectrums_file.open("spectrums.out");
  theoretical_file << "> theoretical" << endl;
  observed_file << "> observed" << endl;
  spectrums_file << "> spectrums" << endl;
  int max_mz = (int)max_mz_float;
  int min_mz = (int)min_mz_float;
  //int max_mz = 1300;
  map<PEAK_T*, string> peak_colors;
  carp(CARP_DEBUG, "min mz: %d, max mz: %d\n", max_mz);
  FLOAT_T average = 0;
  for (int peak_index = 0; peak_index < spectrum->num_peaks; ++peak_index) {
    average = average + get_peak_intensity(find_peak(spectrum->peaks, peak_index));
  }
  average = average / spectrum->num_peaks;
  //cout << "AVERAGE " << average << endl;
  // make spectacle file for observed peaks
  for (int peak_index = 0; peak_index < spectrum->num_peaks; ++peak_index) {
    FLOAT_T location = get_peak_location(find_peak(spectrum->peaks,peak_index));
    FLOAT_T intensity = get_peak_intensity(find_peak(spectrum->peaks, peak_index)); 
    if (location > min_mz && location < max_mz) {
    if (normalize) {
      peak_colors[find_peak(spectrum->peaks, peak_index)] = "blue";
      //observed_file << location<< "\t" << pow(intensity * average * normalize, 0.2) << "\tnolabel\tred" << endl;
      //spectrums_file << location<< "\t" << pow(intensity * average * normalize, 0.2) << "\tnolabel\tblue" << endl;
      //spectrums_file << location<< "\t" << pow(intensity * average * normalize, 0.25) << "\tnolabel" << endl;
    } else {
      observed_file << location<< "\t" << intensity << "\tnolabel\tred" << endl;
      spectrums_file << location<< "\t" << intensity  << "\tnolabel\tblue" << endl;
      //spectrums_file << location<< "\t" << intensity  << "\tnolabel" << endl;
    }
    }
  }
  observed_file.close();
  // make spectacle file for theoretical peaks
  FLOAT_T* index = theoretical;
  int i = 0;
  int match_count = 0;
  int mismatch_count = 0;
  int b_match_count = 0;
  int y_mismatch_count = 0;
  while (i <= max_mz)  {
    if (((*index > 1 && !noflanks) || *index > 26) && i >= min_mz) {
        theoretical_file << i << "\t" << *index << "\tnolabel\tred" << endl;
      PEAK_T* peak = get_nearest_peak(spectrum, i, 1);
      if (peak != NULL) {
	++match_count;
	peak_colors[peak] = "green";
        //if (*index == 50) ++b_match_count;
	spectrums_file << i << "\t" << 40 - *index << "\tnolabel\tgreen" << endl;
	//spectrums_file << get_peak_location(peak) << "\t" << pow (get_peak_intensity(peak) * average * normalize, 0.2) << "\tnolabel\tgreen" << endl;
      } else {
        //if (*index == 51) ++y_mismatch_count;
	++mismatch_count;
        spectrums_file << i << "\t" << 40 - *index << "\tnolabel\tred" << endl;
      }
    }
    ++i;
    ++index;
  }
  FLOAT_T location;
  FLOAT_T intensity;
  for (map<PEAK_T*, string>::iterator it = peak_colors.begin(); it != peak_colors.end(); ++it) {
    location = get_peak_location(it->first);
    intensity = get_peak_intensity(it->first);
    spectrums_file << location << "\t" << pow(intensity * average * normalize, 0.2) << "\tnolabel\t" << it->second << endl;
  }

  cout << "b matches " << b_match_count << " y mismatches " << y_mismatch_count << endl;
  cout << "match: " << match_count << " mismatch: " << mismatch_count << endl;
  theoretical_file.close();
  spectrums_file.close();
}




/*****************************************************/
/* (slightly) modified functions from scorer.c below */
/*****************************************************/
/*
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
  //print_spectrums(theoretical, spectrum, 1200, 1);

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
    ion->calculate_mass();
    intensity_array_idx = (int)(ion->get_mz() / bin_width + 0.5);
    ion_type = ion->type();
    ion_charge = ion->charge();
    //cout << "m/z: " << ion->get_mz() << " charge: " << ion->charge() << endl;
    // skip ions that are located beyond max mz limit
    if(intensity_array_idx >= scorer->sp_max_mz){
      continue;
    }

  //cout << *ion << endl;
  //if (ion->type() == B_ION) 
    //cout << "bion" << endl; else cout << "yion" << endl;

    // is it B, Y ion?

    // neutral loss peak?
    // Add peaks of intensity 50.0 for B, Y type ions. 
    // In addition, add peaks of intensity of 25.0 to +/- 1 m/z flanking each B, Y ion.
    // Skip ions that are located beyond max mz limit
    if((intensity_array_idx)< scorer->sp_max_mz){
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
*/
