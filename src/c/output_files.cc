extern "C" {
#include "carp.h"
#include "parameter.h"
#include "spectrum_collection.h"
#include "match_collection.h"
}
#include "output_files.h"

// C interface
extern "C" void* construct_output_files(int num_decoys, int num_proteins) {
  return static_cast<void*>(new OutputFiles(num_decoys, num_proteins));
}

extern "C" void close_output_files(void* output_files_param, int num_spectra) {
  OutputFiles* output_files = static_cast<OutputFiles*>(output_files_param);
  output_files->Close(num_spectra);
  delete output_files;
}


void OutputFiles::Close(int num_spectra) {
  // fix headers in csm files
  for(int file_idx=0; file_idx < num_decoys_ + 1; file_idx++){
    carp(CARP_DEBUG, "Changing csm header to have %i spectrum searches",
	 num_spectra);
    serialize_total_number_of_spectra(num_spectra,
				      psm_file_array_[file_idx]);
  }

  // No evidence of closing or flushing files in match_search.c.
  // Not sure what's up. -- Benjamin
}

void OutputFiles::PrintMatches(MATCH_COLLECTION_T* match_collection,
			       SPECTRUM_T* spectrum, bool is_decoy,
			       int psm_file_index,
			       bool send_to_sqt_tab) {
  print_matches(match_collection,
		spectrum,
		is_decoy,
		psm_file_array_[psm_file_index],
		sqt_file_,
		send_to_sqt_tab ? decoy_sqt_file_ : NULL,
		tab_file_,
		send_to_sqt_tab ? decoy_tab_file_ : NULL);
}

/**
 * \brief A private function for crux-search-for-matches to prepare
 * binary psm and text sqt files.
 *
 * Reads the --overwrite and --output-mode values from
 * parameter.c. Opens psm file(s) if requested, setting a given
 * pointer to the array of filehandles.  Opens sqt file(s) if
 * requested, setting the given pointers to each file handle.  If
 * binary files not requested, creates an array of NULL pointers.  If
 * sqt files not requested, sets given pointers to NULL. 
 *
 * \returns void.  Sets given arguments to newly created filehandles.
 */
void OutputFiles::Open(int num_proteins) {
  char* match_output_folder = get_string_parameter("match-output-folder");
  MATCH_SEARCH_OUTPUT_MODE_T output_type = get_output_type_parameter(
                                                    "output-mode");
  BOOLEAN_T overwrite = get_boolean_parameter("overwrite");
  carp(CARP_DEBUG, "The output type is %d (binary, sqt, tab, all)" \
       " and overwrite is '%d'", (int)output_type, (int)overwrite);


  // create binary psm files (allocate memory, even if not used)
  psm_file_array_ = create_psm_files();

  if(output_type == SQT_OUTPUT || output_type == ALL_OUTPUT){

    //create sqt file handles
    carp(CARP_DEBUG, "Opening sqt files");
    char* sqt_filename = get_string_parameter_pointer("sqt-output-file");
    sqt_file_ = create_file_in_path(sqt_filename, 
				    match_output_folder, 
				    overwrite);
    char* decoy_sqt_filename = get_string_parameter_pointer(
                                                    "decoy-sqt-output-file");
    if( get_int_parameter("number-decoy-set") > 0 ){
      decoy_sqt_file_ = create_file_in_path(decoy_sqt_filename,
                                            match_output_folder,
                                            overwrite);
    }

    if(sqt_file_ == NULL || decoy_sqt_file_ == NULL){
      carp(CARP_DEBUG, "sqt file or decoy is null");
    }
  }

  if(output_type == TAB_OUTPUT || output_type == ALL_OUTPUT) {

    //create sqt file handles
    carp(CARP_DEBUG, "Opening tab delimited files");
    char* tab_filename = get_string_parameter_pointer("tab-output-file");
    tab_file_ = create_file_in_path(tab_filename, 
                                    match_output_folder, 
                                    overwrite);
    char* decoy_tab_filename = get_string_parameter_pointer(
                                                    "decoy-tab-output-file");
    if( get_int_parameter("number-decoy-set") > 0 ){
      decoy_tab_file_ = create_file_in_path(decoy_tab_filename,
                                            match_output_folder,
                                            overwrite);
    }

    if(tab_file_ == NULL || decoy_tab_file_ == NULL){
      carp(CARP_DEBUG, "tab file or decoy tab file is null");
    }

  }

  free(match_output_folder);
  carp(CARP_DEBUG, "Finished opening output files");

  //print headers
  serialize_headers(psm_file_array_);
  print_sqt_header(sqt_file_, "target", num_proteins, FALSE);// !analyze-matches
  print_sqt_header(decoy_sqt_file_, "decoy", num_proteins, FALSE);
  print_tab_header(tab_file_);
  print_tab_header(decoy_tab_file_);
}
