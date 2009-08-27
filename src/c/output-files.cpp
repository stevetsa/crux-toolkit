/**
 * \file output-files.cpp
 */
/*
 * FILE: output-files.h
 * AUTHOR: Barbara Frewen
 * CREATE DATE: Aug 24, 2009
 * PROJECT: crux
 * DESCRIPTION: A class description for handling all the various
 * output files, excluding parameter and log files.  The filenames,
 * locations and overwrite status would be taken from parameter.c.
 */

#include "output-files.h"
using namespace std;

/**
 * Default constructor for OutputFiles.  Opens all of the needed
 * files, naming them based on the values of the parameters output-dir
 * and fileroot and on the name given (search, percolator, etc.).
 * Requires that the output directory already exist. 
 */
OutputFiles::OutputFiles(const char* program_name){

  psm_file_array_ = NULL;
  tab_file_array_ = NULL;
  sqt_file_array_ = NULL;

  // parameters for all three file types
  BOOLEAN_T overwrite = get_boolean_parameter("overwrite");
  const char* output_directory = get_string_parameter_pointer("output-dir");
  const char* fileroot = get_string_parameter_pointer("fileroot");
  if( strcmp(fileroot, "__NULL_STR") == 0 ){
    fileroot = NULL;
  }
  int num_decoy_files = get_int_parameter("num-decoy-files");
  num_files_ = num_decoy_files + 1; // plus target file

  carp(CARP_DEBUG, 
       "OutputFiles is opening %d files (%d decoys) in '%s' with root '%s'."
       " Overwrite: %d.", 
       num_files_, num_decoy_files, output_directory, fileroot, overwrite);

  // all operations create tab files
  createFiles(&tab_file_array_, 
              output_directory, 
              fileroot, 
              program_name, 
              "txt", 
              overwrite); 

  // search operations create .csm files
  if( strcmp(program_name, "search") == 0 ||
      strcmp(program_name, "sequest") == 0 ){
    createFiles(&psm_file_array_, 
                 output_directory, 
                 fileroot, 
                 program_name, 
                 "csm", 
                 overwrite);
  }

  // only sequest creates sqt files
  if( strcmp(program_name, "sequest") == 0 ){
    createFiles(&sqt_file_array_, 
                 output_directory, 
                 fileroot, 
                 program_name, 
                 "sqt", 
                 overwrite);
  }
}

OutputFiles::~OutputFiles(){
  delete [] psm_file_array_;
  delete [] tab_file_array_;
  delete [] sqt_file_array_;
}

/**
 * A private function for generating target and decoy files named
 * according to the given arguments.
 *
 * New files are returned via the file_array_ptr argument.  When
 * num_files > 1, exactly one target file is created and the remaining
 * are decoys.  Files are named 
 * "output-dir/fileroot.command_name.target|decoy[n].extension".
 * Requires that the output-dir already exist and have write
 * permissions. 
 * \returns TRUE if num_files new files are created, else FALSE.
 */
BOOLEAN_T OutputFiles::createFiles(FILE*** file_array_ptr,
                                   const char* output_dir,
                                   const char* fileroot,
                                   const char* command_name,
                                   const char* extension,
                                   BOOLEAN_T overwrite){

  if( num_files_ == 0 ){
    return FALSE;
  }
  
  // allocate array
  *file_array_ptr = (FILE**)mycalloc(num_files_, sizeof(FILE*));

  // determine the name for each of the files: target, decoy-1, etc.
  string* target_decoy = new string[num_files_];
  target_decoy[0] = "target";
  if( num_files_ == 2 ){
    target_decoy[1] = "decoy";
  }else{
    for(int file_idx = 1; file_idx < num_files_; file_idx++){
      ostringstream name_builder;
      name_builder << "decoy-" << file_idx;
      target_decoy[file_idx] = name_builder.str();
    }
  }
  
  // create each file
  for(int file_idx = 0; file_idx < num_files_; file_idx++ ){
    // concatinate the pieces of the name
    ostringstream name_builder;
    if( fileroot ){
      name_builder << fileroot << "." ;
    }
    name_builder << command_name << "."
                 << target_decoy[file_idx] << "." 
                 << extension;
    string filename = name_builder.str();

    // open the file (it checks for success)
    (*file_array_ptr)[file_idx] = create_file_in_path(filename.c_str(),
                                                    output_dir,
                                                    overwrite);
  }// next file

  delete [] target_decoy;

  return TRUE;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
























