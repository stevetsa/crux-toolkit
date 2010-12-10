#include "spectral-counts.h"

using namespace std;

void getDirPath(char* path, char** dir);
set<MATCH_T*> filterMatches(MATCH_COLLECTION_ITERATOR_T* match_collection_it);
map<PEPTIDE_T*, FLOAT_T> getPeptideScores();

int spectral_counts_main(int argc, char** argv){
  const char* option_list[] = {
    "threshold",
    "input-ms2",
    "fileroot",
    "output-dir",
    "overwrite",
    "unique-mapping",
    "input-bullseye",
    "quant-level",
    "measure",
    "average",
    "version",
    "verbosity"
  };
  const char* argument_list[] = {
    "input-PSM",
    "database"
  };

  int num_options = sizeof(option_list) / sizeof(char*);
  int num_arguments = sizeof(argument_list) / sizeof(char*);


  initialize_run(SPECTRAL_COUNTS_COMMAND, argument_list, num_arguments,
		 option_list, num_options, argc, argv);

  /* Create Match Collection from input-PSM */
  char * psm_file = get_string_parameter("input-PSM");
  char * database = get_string_parameter("database");
  char * dir_path = NULL;

  getDirPath(psm_file, &dir_path);
  
  int decoy_count = 0;

  MATCH_COLLECTION_ITERATOR_T* match_collection_it 
    = new_match_collection_iterator(dir_path, database, &decoy_count);
  
  set<MATCH_T*> matches = filterMatches(match_collection_it);
  carp(CARP_INFO, "Number of matches passed the threshold %i", matches.size());
  
  map<PEPTIDE_T*, FLOAT_T> peptideToScore = getPeptideScores();
  
  


  free(dir_path);
  return 1;
}


map<PEPTIDE_T*, FLOAT_T> getPeptideScores(){
  map<PEPTIDE_T*, FLOAT_T> peptide_scores;
  char * measure = get_string_parameter("measure"); //TODO cehck this value
  char * ms2 = get_string_parameter("input-MS2");

  if (measure == NULL || strcmp(measure, "SIN")==0){
    map<pair<int,int>, Spectrum*> spectras;
    SpectrumCollection* spectrumCollection = new SpectrumCollection(ms2);
    if (!spectrumCollection->parse()){
      carp(CARP_FATAL, "Failed to parse ms2 file: %s", ms2);
    } 
    FilteredSpectrumChargeIterator* spectrum_iterator = 
      new FilteredSpectrumChargeIterator(spectrumCollection);

    while (spectrum_iterator->hasNext()){
      int charge = 0;
      Spectrum* spectrum = spectrum_iterator->next(&charge);
      spectras.insert(make_pair(make_pair(spectrum->getFirstScan(), 
	spectrum->getLastScan()), spectrum));
    }
    
  }
  return peptide_scores;
}




/*
 * Helper function to find the path of directory
 * for a given path to a file
 * TODO: write unit tests
 */
void getDirPath(char* path, char** dir){
  int pos = strlen(path)-1;
  while (pos >0 && *(path+pos) != '/'){pos--;}
  if (pos == 0){ // no backslashes exist
    *dir = (char*) malloc(sizeof(char)*2);
    strcpy(*dir, ".");
  } else {
    *dir = (char*) malloc((sizeof(char)*pos)+1);
    strncpy(*dir, path, pos);
    (*dir)[pos] = '\0';
  }
}


set<MATCH_T*> filterMatches(MATCH_COLLECTION_ITERATOR_T* match_collection_it){
  set<MATCH_T*> matches;
  MATCH_ITERATOR_T* match_iterator = NULL;
  MATCH_COLLECTION_T* match_collection = NULL;
  FLOAT_T threshold = get_double_parameter("threshold");
  BOOLEAN_T qualify = FALSE;
  while (match_collection_iterator_has_next(match_collection_it)){

    
    match_collection = match_collection_iterator_next(match_collection_it);
    match_iterator = new_match_iterator(match_collection, XCORR, TRUE);

    while(match_iterator_has_next(match_iterator)){
      MATCH_T* match = match_iterator_next(match_iterator);
      qualify = FALSE;
      if (get_match_rank(match, XCORR) != 1){
	continue;
      }
      // find a qvalue score lower than threshold
      if (get_match_score(match, PERCOLATOR_QVALUE) != FLT_MIN &&
	  get_match_score(match, PERCOLATOR_QVALUE) <= threshold)  {
	qualify = TRUE;
      } else if (get_match_score(match, QRANKER_QVALUE) != FLT_MIN ||
	       get_match_score(match, QRANKER_QVALUE) <= threshold)  {
	qualify = TRUE;
      } else if (get_match_score(match, DECOY_XCORR_QVALUE) != FLT_MIN ||
		 get_match_score(match, DECOY_XCORR_QVALUE) <= threshold)  {
	qualify = TRUE;
      } 
      
      if (qualify == TRUE){
	matches.insert(match);
      }
    }
  }
  return matches;
}
