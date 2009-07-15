/*************************************************************************//**
 * \file spectrum.cpp
 * AUTHOR: Chris Park
 * CREATE DATE:  June 22 2006
 * DESCRIPTION: code to support working with spectra
 * REVISION: $Revision: 1.71 $
 ****************************************************************************/
#include "spectrum.h"
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#include "objects.h"

#include "peak.h"
#include "utils.h"
#include "mass.h"
#include "parameter.h"

#include "carp.h"

/**
 * \define constants
 */
#define MZ_TO_PEAK_ARRAY_RESOLUTION 5 // i.e. 0.2 m/z unit
#define MAX_PEAK_MZ 5000
#define MAX_CHARGE 4


/***************************************************
 *SPECTRUM_T: PROTECTED MEMBER FUNCTIONS
 ***************************************************/


/***************************************************
 *SPECTRUM_T: PUBLIC MEMBER FUNCTIONS
 ***************************************************/
void SPECTRUM_T::init() {
  //carp(CARP_INFO,"SPECTRUM_T::init(): start.");
  first_scan = 0;
  last_scan = 0;
  id = 0;
  spectrum_type = (SPECTRUM_TYPE_T)0;
  precursor_mz = 0;
  num_possible_z = 0;
  peaks.clear();
  min_peak_mz = 0;
  max_peak_mz = 0;
  total_energy = 0;
  filename = NULL;
  i_lines[0] = '\0';
  d_lines[0] = '\0';
  has_peaks = FALSE;
  sorted_by_mz = FALSE;
  sorted_by_intensity = FALSE;
  has_mz_peak_array = FALSE;
  mz_peak_array = NULL;


  int line_idx;
  possible_z = (int*)mymalloc(sizeof(int) * MAX_CHARGE);
  
  // initialize D lines
  for(line_idx = 0; line_idx < MAX_D_LINES; ++line_idx){
    d_lines[line_idx] = NULL;
  }

  // initialize I lines
  for(line_idx = 0; line_idx < MAX_I_LINES; ++line_idx){
    i_lines[line_idx] = NULL;
  }  
  has_peaks = FALSE;
  //carp(CARP_INFO,"SPECTRUM_T::init() end");

}
/**
 * \creates An (empty) spectrum object.
 */
SPECTRUM_T::SPECTRUM_T() {
  init();
}



/**
 * \creates A new spectrum object, populated with the user specified parameters.
 */ 
SPECTRUM_T::SPECTRUM_T(
  int               first_scan,         ///< The number of the first scan -in
  int               last_scan,          ///< The number of the last scan -in
  SPECTRUM_TYPE_T   spectrum_type,      ///< The type of spectrum. -in
  FLOAT_T             precursor_mz,       ///< The m/z of the precursor (for MS-MS spectra) -in
  int*              possible_z,         ///< The possible charge states of this spectrum  -in
  int               num_possible_z,     ///< The number of possible charge states of this spectrum  -in
  char*             filename)           ///< Optional filename -in
{
  init();
  
  this->first_scan = first_scan;
  this->last_scan = last_scan;
  this->spectrum_type = spectrum_type;
  this->precursor_mz = precursor_mz;
  this->sorted_by_mz = FALSE;
  this->has_mz_peak_array = FALSE;
  this->sorted_by_intensity = FALSE;
  this->mz_peak_array = NULL;
  set_new_possible_z(possible_z, num_possible_z);
  set_new_filename(filename);
}
/*****************************************
 *Creates a spectrum_t object from a MSToolkit::Spectrum object
 *****************************************/
SPECTRUM_T::SPECTRUM_T(MSToolkit::Spectrum& mst_spectrum) {
  init();
  first_scan = mst_spectrum.getScanNumber();
  last_scan = mst_spectrum.getScanNumber();


  spectrum_type = MS2; //assume MS2
  switch (mst_spectrum.getFileType()) {
  case MSToolkit::MS1:
    spectrum_type = MS1;
    break;
  case MSToolkit::MS2:
    spectrum_type = MS2;
    break;
  case MSToolkit::MS3:
    spectrum_type = MS3;
    break;
  case MSToolkit::ZS:
  case MSToolkit::UZS:
  case MSToolkit::IonSpec:
  case MSToolkit::SRM:
  case MSToolkit::Unspecified:
  default:
    carp(CARP_DETAILED_DEBUG,"Unsuppported spectra type!:%d",mst_spectrum.getFileType());
    carp(CARP_DETAILED_DEBUG,"Assuming ms2");
  }
  
  precursor_mz = mst_spectrum.getMZ();

  for (int z_idx=0;z_idx < mst_spectrum.sizeZ();z_idx++) {
    add_possible_z(mst_spectrum.atZ(z_idx).z);
  }
  
  for (int peak_idx=0;peak_idx < mst_spectrum.size();peak_idx++) {
    add_peak(mst_spectrum[peak_idx].intensity,
	     mst_spectrum[peak_idx].mz);
  }
}

/**
 * Destructor for a spectrum object.
 */
SPECTRUM_T::~SPECTRUM_T() {
  int line_idx;
  
  // only non post_process spectrum has these features to free
  if(has_peaks){
    free(possible_z);
    peaks.clear();
    free(filename);
    
    // free D lines
    for(line_idx = 0; line_idx < MAX_D_LINES; ++line_idx){
      if(d_lines[line_idx] != NULL){
        free(d_lines[line_idx]);
      }
      else{
        break;
      }
    }
    
    // free I lines
    for(line_idx = 0; line_idx < MAX_I_LINES; ++line_idx){
      if(i_lines[line_idx] != NULL){
        free(i_lines[line_idx]);
      }
      else{
        break;
      }
    }    
  }
  
  if(has_mz_peak_array){
    free(mz_peak_array);
  }
}

/**
 * Prints a spectrum object to file.
 */
void SPECTRUM_T::print(
  FILE* file ///< output file to print at -out
  )
{
  int num_z_index = 0;
  int num_d_index = 0;
  int num_i_index = 0;
  unsigned int num_peak_index = 0;

  fprintf(file, "S\t%06d\t%06d\t%.2f\n", 
         first_scan,
         last_scan,
         precursor_mz);

  // print 'I' line
  for(; num_i_index < MAX_I_LINES; ++num_i_index){
    if(i_lines[num_i_index] == NULL){
      break;
    }

    fprintf(file, "%s", i_lines[num_i_index]);
  }
  
  // print 'Z', 'D' line
  for(; num_z_index < num_possible_z; ++num_z_index){

    // print 'Z' line
    fprintf(file, "Z\t%d\t%.2f\n", possible_z[num_z_index],
            get_singly_charged_mass(possible_z[num_z_index]));

    // are there any 'D' lines to print?
    if(num_d_index < MAX_D_LINES){
      if(d_lines[num_d_index] != NULL){
        fprintf(file, "%s", d_lines[num_d_index]);
      }
      ++num_d_index;
    }
  }
  
  // print peaks
  for(; num_peak_index < peaks.size(); ++num_peak_index){
    fprintf(file, "%.2f %.13f\n",
	    peaks[num_peak_index].get_location(),
	    peaks[num_peak_index].get_intensity());
  }
}

void SPECTRUM_T::print_sqt(
  FILE* file,           ///< output file to print at -out
  int num_matches,      ///< number of peptides compared to this spec -in
  int charge            ///< charge used for the search -in
  ){

  int precision = get_int_parameter("precision");
  char format[64];
  sprintf(format, 
          "S\t%%d\t%%d\t%%d\t%%.%if\t%%s\t%%.%if\t%%.%if\t%%.%if\t%%d\n", 
          precision, precision, precision, precision);
  //<first scan><last scan><charge><precursor m/z><# sequence match>
  //fprintf(file, "S\t%d\t%d\t%d\t%.2f\t%s\t%.2f\t%.2f\t%.2f\t%d\n", 
  fprintf(file, format,
          get_first_scan(),
          get_last_scan(),
          charge, 
          0.0, // FIXME dummy <process time>
          "server", // FIXME dummy <server>
          // get_spectrum_precursor_mz(spectrum), // TODO this should be M+H+
          get_neutral_mass(charge), //this is used in search
          0.0, // FIXME dummy <total intensity>
          0.0, // FIXME dummy <lowest sp>
          num_matches);
  
}


/**
 * Prints a spectrum object to STDOUT.
 */
void SPECTRUM_T::print_stdout(
  )
{
  print(stdout);
}

/**
 * Copies spectrum object src to dest.
 * must pass in a memory allocated SPECTRUM_T* dest
 * doesn't copy the sum array related fields
 */
void SPECTRUM_T::copy(
  SPECTRUM_T* src, ///< the source spectrum -in
  SPECTRUM_T* dest ///< the destination spectrum -out
  )
{
  int num_peak_index = 0;
  int* possible_z;
  char* new_filename;
  int line_idx;
  
  // copy each varible
  dest -> set_first_scan(src -> get_first_scan());
  dest -> set_last_scan(src -> get_last_scan());
  dest -> set_id(src -> get_id());
  dest -> set_spectrum_type(src -> get_spectrum_type());
  dest -> set_precursor_mz(src -> get_precursor_mz());
  
  // copy possible_z
  possible_z = src -> get_possible_z();
  dest -> set_possible_z(possible_z, 
                          src -> get_num_possible_z());
  free(possible_z);
  
  // copy filename
  new_filename = src -> get_filename();
  dest -> set_filename(new_filename);
  free(new_filename);
  
  // copy 'D', 'I' lines
  for(line_idx = 0; line_idx < MAX_D_LINES; ++line_idx){
    if(src->d_lines[line_idx] != NULL){
       dest->d_lines[line_idx] = my_copy_string(src->d_lines[line_idx]);
    }
    else{
      break;
    }
  }
  
  // copy 'D', 'I' lines
  for(line_idx = 0; line_idx < MAX_I_LINES; ++line_idx){
    if(src->i_lines[line_idx] != NULL){
       dest->i_lines[line_idx] = my_copy_string(src->i_lines[line_idx]);
    }
    else{
      break;
    }
  }

  // copy each peak
  for(; num_peak_index < src -> get_num_peaks(); ++num_peak_index){
    dest -> add_peak(src -> peaks[num_peak_index].get_intensity(),
		     src -> peaks[num_peak_index].get_location());
  }
}

/**
 * Parses a spectrum from file.
 * \returns TRUE if success. FALSE is failure.
 * 'I'
 * Skips Header line "H"
 * FIXME if need to read 'H', header line, does not parse ID
 */
BOOLEAN_T SPECTRUM_T::parse_spectrum_file(
  FILE* file, ///< the input file stream -in
  char* filename ///< filename of the spectrum, should not free -in
  )
{
  long file_index = ftell(file); // stores the location of the current working line in the file
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  FLOAT_T location_mz;
  FLOAT_T intensity;
  BOOLEAN_T record_S = FALSE; // check's if it read S line
  BOOLEAN_T record_Z = FALSE; // check's if it read Z line
  BOOLEAN_T start_add_peaks = FALSE; // check's if it started reading peaks
  BOOLEAN_T file_format = FALSE; // is the file format correct so far
  
  FLOAT_T test_float;
  char test_char;
  
  while( (line_length = getline(&new_line, &buf_length, file)) != -1){
    // skip header line
    // if(new_line[0] == 'H'){
    //  file_index = ftell(file);
    //  continue;
    // }
    // checks if 'S' is not the first line
    if((!record_S || (record_S && start_add_peaks)) && 
            (new_line[0] == 'Z' ||  
             new_line[0] == 'I' ||
             new_line[0] == 'D' )){
      file_format = FALSE;
      carp(
	   CARP_ERROR, 
	   "Incorrect order of line (S,Z, Peaks)\n"
	   "At line: %s", 
	   new_line
	   );
      break; // File format incorrect
    }
    // Reads the 'S' line
    else if(new_line[0] == 'S' && !record_S){
      record_S = TRUE;
      if(!parse_S_line(new_line, buf_length)){
        file_format = FALSE;
        break; // File format incorrect
      }
    }
    // Reads the 'Z' line 
    else if(new_line[0] == 'Z'){
      record_Z = TRUE;
      if(!parse_Z_line(new_line)){
        file_format = FALSE;
        break; // File format incorrect
      }
    }
    
    // Reads the 'D' line 
    else if(new_line[0] == 'D'){
      if(!parse_D_line(new_line)){
        file_format = FALSE;
        break; // File format incorrect
      }
    }
    
    // Reads the 'I' line 
    else if(new_line[0] == 'I'){
      if(!parse_I_line( new_line)){
        file_format = FALSE;
        break; // File format incorrect
      }
    }
    
    // Stops, when encounters the start of next spectrum 'S' line
    else if(new_line[0] == 'S' && start_add_peaks){ // start of next spectrum
      break;
    }
    
    // *****parse peak line******
    else if(new_line[0] != 'Z' &&  
            new_line[0] != 'I' &&
            new_line[0] != 'D' &&
            new_line[0] != '\n')
      {
        // checks if the peaks are in correct order of lines
        if((!record_Z || !record_S)){
          file_format = FALSE;
          carp(
	       CARP_ERROR,
	       "Incorrect order of line (S,Z, Peaks)\n"
	       "At line: %s", 
	       new_line
	       );
          break; // File format incorrect
        }
        // check for peak line format
#ifdef USE_DOUBLES
        else if((sscanf(new_line,"%lf %lf %lf",// test format:peak line has more than 2 fields
                        &test_float, &test_float, &test_float) > 2)||
                (sscanf(new_line,"%lf %lf %c",// test format:peak line has more than 2 fields
                        &test_float, &test_float, &test_char) > 2)||
                (sscanf(new_line,"%lf %lf",// test format:peak line has less than 2 fields
                        &test_float, &test_float) != 2)){
#else
        else if((sscanf(new_line,"%f %f %f",// test format:peak line has more than 2 fields
                        &test_float, &test_float, &test_float) > 2)||
                (sscanf(new_line,"%f %f %c",// test format:peak line has more than 2 fields
                        &test_float, &test_float, &test_char) > 2)||
                (sscanf(new_line,"%f %f",// test format:peak line has less than 2 fields
                        &test_float, &test_float) != 2)){
#endif
          file_format = FALSE;
          carp(
	       CARP_ERROR,
	       "Incorrect peak line\n"
	       "At line: %s", 
	       new_line
	       );
          break; // File format incorrect
        }
        // Reads the 'peak' lines, only if 'Z','S' line has been read
#ifdef USE_DOUBLES
        else if(record_Z && record_S &&
                (sscanf(new_line,"%lf %lf", &location_mz, &intensity) == 2))
#else
        else if(record_Z && record_S &&
                (sscanf(new_line,"%f %f", &location_mz, &intensity) == 2))
#endif
	  {
	    file_format = TRUE;
	    start_add_peaks = TRUE;
	    add_peak(intensity, location_mz);
	  }
	}
	//*************************
	    file_index = ftell(file); // updates the current working line location
      }
  

    // now we have peak information
    has_peaks = TRUE;
    
    // set the file pointer back to the start of the next 's' line
    fseek(file, file_index, SEEK_SET);
    myfree(new_line);
    
    // set filename of empty spectrum
    set_new_filename(filename);
    
    // No more spectrum in .ms file
    if(!record_S && !file_format){
      return FALSE;
    }
    
    // File format incorrect
    if(!file_format){ 
      carp(CARP_ERROR, "incorrect file format\n");
      return FALSE;
    }
    return TRUE;
}


/**
 * Parses the 'S' line of the a spectrum
 * \returns TRUE if success. FALSE is failure.
 * 
 */
BOOLEAN_T SPECTRUM_T::parse_S_line(
   char* line, ///< 'S' line to parse -in
   int buf_length ///< line length -in
   )
{
  char spliced_line[buf_length];
  int line_index = 0;
  int spliced_line_index = 0;
  int first_scan;
  int last_scan;
  FLOAT_T precursor_mz;
  FLOAT_T test_float;
  char test_char;
  
  // deletes empty space & 0
  while((line[line_index] !='\0') && 
	(line[line_index] == 'S' || 
	 line[line_index] == '\t'||
	 line[line_index] == ' ' || 
	 line[line_index] == '0')){
    ++line_index;
  }
  // reads in line value
  while(line[line_index] !='\0' && 
	line[line_index] != ' ' && 
	line[line_index] != '\t'){
    spliced_line[spliced_line_index] =  line[line_index];
    ++spliced_line_index;
     ++line_index;
  }
   spliced_line[spliced_line_index] =  line[line_index];
   ++spliced_line_index;
   ++line_index;
   // deletes empty space & zeros
   while((line[line_index] !='\0') && 
         (line[line_index] == '\t' || 
          line[line_index] == ' ' || 
          line[line_index] == '0')){
     ++line_index;
   }
   // read last scan & precursor m/z
   while(line[line_index] !='\0'){
     spliced_line[spliced_line_index] =  line[line_index];
     ++spliced_line_index;
     ++line_index;
   }
   spliced_line[spliced_line_index] = '\0';
   
   // check if S line is in correct format
#ifdef USE_DOUBLES
   if ( (sscanf(spliced_line,"%lf %lf %lf %lf",// test format:S line has more than 3 fields
                &test_float, &test_float, &test_float, &test_float) > 3) ||
        (sscanf(spliced_line,"%lf %lf %lf %c",// test format:S line has more than 3 fields 
                &test_float, &test_float, &test_float, &test_char) > 3) ||
        (sscanf(spliced_line,"%i %i %lf", // S line is parsed here
		&first_scan, &last_scan, &precursor_mz) != 3)) 
#else
     if ( (sscanf(spliced_line,"%f %f %f %f",// test format:S line has more than 3 fields
		  &test_float, &test_float, &test_float, &test_float) > 3) ||
	  (sscanf(spliced_line,"%f %f %f %c",// test format:S line has more than 3 fields 
		  &test_float, &test_float, &test_float, &test_char) > 3) ||
	  (sscanf(spliced_line,"%i %i %f", // S line is parsed here
		&first_scan, &last_scan, &precursor_mz) != 3)) 
#endif
       {
	 carp(CARP_ERROR,"Failed to parse 'S' line:\n %s",line);
	 return FALSE;
       }
     
   set_first_scan(first_scan);
   set_last_scan(last_scan);
   set_precursor_mz(precursor_mz);
   
   return TRUE;
}

/**
 * Parses the 'Z' line of the a spectrum
 * \returns TRUE if success. FALSE is failure.
 * 
 */
BOOLEAN_T SPECTRUM_T::parse_Z_line(
   char* line  ///< 'Z' line to parse -in
   )
{
  int tokens;
  char line_name;
  int charge;
  FLOAT_T m_h_plus;
  FLOAT_T test_float;
  char test_char;
  
  // check if Z line is in correct format
#ifdef USE_DOUBLES
  if( ((tokens =  // test format: Z line has less than 3 fields
	sscanf(line, "%c %lf %lf", &test_char, &test_float, &test_float)) < 3) ||
      ((tokens =   // test format: Z line has more than 3 fields
	sscanf(line, "%c %lf %lf %lf", &test_char, &test_float, &test_float, &test_float)) >  3) ||
      ((tokens =  // test format: Z line has more than 3 fields
	sscanf(line, "%c %lf %lf %c", &test_char, &test_float, &test_float, &test_char)) >  3) ||
      (tokens = // Z line is parsed here
       sscanf(line, "%c %d %lf", &line_name, &charge, &m_h_plus)) != 3)
#else
    if( ((tokens =  // test format: Z line has less than 3 fields
         sscanf(line, "%c %f %f", &test_char, &test_float, &test_float)) < 3) ||
       ((tokens =   // test format: Z line has more than 3 fields
         sscanf(line, "%c %f %f %f", &test_char, &test_float, &test_float, &test_float)) >  3) ||
       ((tokens =  // test format: Z line has more than 3 fields
         sscanf(line, "%c %f %f %c", &test_char, &test_float, &test_float, &test_char)) >  3) ||
       (tokens = // Z line is parsed here
        sscanf(line, "%c %d %f", &line_name, &charge, &m_h_plus)) != 3)
#endif
      {
	carp(CARP_ERROR,"Failed to parse 'Z' line:\n %s",line);
	return FALSE;
      }  
  
  return add_possible_z(charge);
}

/**
 * Adds a possible charge(z) to the spectrum.
 * Must not exceed the MAX_CARGE capacity
 */
BOOLEAN_T SPECTRUM_T::add_possible_z(
   int charge  ///< charge to add
   )
 {
   int* possible_charge = (int *)mymalloc(sizeof(int));
   *possible_charge = charge;
   if(num_possible_z < MAX_CHARGE){ // change to dynamic sometime...
     possible_z[num_possible_z] = *possible_charge; 
     ++num_possible_z;
     free(possible_charge);
     return TRUE;
   }
   free(possible_charge);
   return FALSE;
 }


/**
 * FIXME currently does not parse D line, just copies the entire line
 * Parses the 'D' line of the a spectrum
 * \returns TRUE if success. FALSE is failure.
 */
BOOLEAN_T SPECTRUM_T::parse_D_line(
   char* line  ///< 'D' line to parse -in
   )
{
  int line_idx;
  int length = strlen(line)+1;
  char* d_line = (char*)mycalloc(length, sizeof(char));
  
  strncpy(d_line, line, length-3);
  d_line[length-2] = '\0';
  d_line[length-3] = '\n';
  
  // find empty spot D lines
  for(line_idx = 0; line_idx < MAX_D_LINES; ++line_idx){
    // check for empty space
    if(d_lines[line_idx] == NULL){
      d_lines[line_idx] = d_line;
      break;
    }
  }

  // check if added new d line to spectrum
  if(line_idx == MAX_D_LINES){
    free(d_line);
    carp(CARP_WARNING, "no more space for additional D lines, max: %d", MAX_D_LINES);
  }
  
  return TRUE;
}

/**
 * FIXME currently does not parse I line, just copies the entire line
 * Parses the 'I' line of the a spectrum
 * \returns TRUE if success. FALSE is failure.
 */
BOOLEAN_T SPECTRUM_T::parse_I_line(
   char* line  ///< 'I' line to parse -in
   )
 {
   int line_idx;
   int length = strlen(line)+1;
   char* i_line = (char*)mycalloc(length, sizeof(char));
   
   strncpy(i_line, line, length-3);
   i_line[length-2] = '\0';
   i_line[length-3] = '\n';
   
   // find empty spot I lines
   for(line_idx = 0; line_idx < MAX_I_LINES; ++line_idx){
     // check for empty space
    if(i_lines[line_idx] == NULL){
      i_lines[line_idx] = i_line;
      break;
    }
  }

  // check if added new i line to spectrum
  if(line_idx == MAX_I_LINES){
    free(i_line);
    carp(CARP_WARNING, "no more space for additional I lines, max: %d", MAX_I_LINES);
  }

  return TRUE;
}

/**
 * Adds a peak to the spectrum given a intensity and location
 * calls update_spectrum_fields to update num_peaks, min_peak ...
 */
BOOLEAN_T SPECTRUM_T::add_peak(
  FLOAT_T intensity, ///< the intensity of peak to add -in
  FLOAT_T location_mz ///< the location of peak to add -in
  )
{
  PEAK_T peak(intensity, location_mz);
  peaks.push_back(peak);
  update_spectrum_fields(intensity, location_mz);
  return TRUE;
}

void SPECTRUM_T::populate_mz_peak_array()
{
  if (has_mz_peak_array == TRUE){
    return;
  }
  
  int array_length = MZ_TO_PEAK_ARRAY_RESOLUTION * MAX_PEAK_MZ;
  PEAK_T** mz_peak_array = (PEAK_T**)mymalloc(array_length * sizeof(PEAK_T*));
  int peak_idx;
  for (peak_idx = 0; peak_idx < array_length; peak_idx++){
    mz_peak_array[peak_idx] = NULL;
  }

  for (PEAK_ITERATOR_T peak_iterator = begin();
       peak_iterator != end();
       ++peak_iterator) {
    FLOAT_T peak_mz = peak_iterator -> get_location();
    int mz_idx = (int) (peak_mz * MZ_TO_PEAK_ARRAY_RESOLUTION);
    if (mz_peak_array[mz_idx] != NULL){
      carp(CARP_INFO, "Peak collision at mz %.3f = %i", peak_mz, mz_idx);
      if(mz_peak_array[mz_idx] -> get_intensity()< peak_iterator -> get_intensity()){
        mz_peak_array[mz_idx] = &(*peak_iterator);
      }
    } else {
      mz_peak_array[mz_idx] = &(*peak_iterator); 
    }
  }
  this->mz_peak_array = mz_peak_array;
  this->has_mz_peak_array = TRUE;
}

/**
 * \returns The closest intensity within 'max' of 'mz' in 'spectrum'
 * NULL if no peak.
 * This should lazily create the data structures within the
 * spectrum object that it needs.
 * TODO: reimplement with faster peak lookup
 */
PEAK_T* SPECTRUM_T::get_nearest_peak(
  FLOAT_T mz, ///< the mz of the peak around which to sum intensities -in
  FLOAT_T max ///< the maximum distance to get intensity -in
  )
{
  populate_mz_peak_array(); // for rapid peak lookup by mz

  FLOAT_T min_distance = BILLION;
  int min_mz_idx = (int)((mz - max) * MZ_TO_PEAK_ARRAY_RESOLUTION + 0.5);
  min_mz_idx = min_mz_idx < 0 ? 0 : min_mz_idx;
  int max_mz_idx = (int)((mz + max) * MZ_TO_PEAK_ARRAY_RESOLUTION + 0.5);
  int absolute_max_mz_idx = MAX_PEAK_MZ * MZ_TO_PEAK_ARRAY_RESOLUTION - 1;
  max_mz_idx = max_mz_idx > absolute_max_mz_idx 
    ? absolute_max_mz_idx : max_mz_idx;
  PEAK_T* peak = NULL;
  PEAK_T* nearest_peak = NULL;
  int peak_idx;
  for (peak_idx=min_mz_idx; peak_idx < max_mz_idx + 1; peak_idx++){
    if ((peak = mz_peak_array[peak_idx]) == NULL){
      continue;
    }
    FLOAT_T peak_mz = peak -> get_location();
    //FLOAT_T distance = abs((FLOAT_T)mz - (FLOAT_T)peak_mz);
    FLOAT_T distance = mz - peak_mz;
    if (distance < 0) distance = -distance;
    if (distance > max){
      continue;
    }
    if (distance < min_distance){
      nearest_peak = peak;
      min_distance = distance;
    }
  }
  return nearest_peak;
}

/**
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T SPECTRUM_T::parse_spectrum(
  char*      filename ///< the file to parse -in
  )
{
  FILE* file;
  if ((file = fopen(filename,"r")) == NULL) {
    carp(CARP_FATAL,"File %s could not be opened",filename);
    return (FALSE);  // exit(1);
  }
  // might check if spectrum is NULL?? Close file??
  if(parse_spectrum_file(file, filename)){
    fclose(file);
    return TRUE;
  }
  fclose(file);
  return FALSE;
}

/**
 * updates num_peaks, min_peak_mz, max_peak_mz, total_energy fields in spectrum
 */
void SPECTRUM_T::update_spectrum_fields(
  FLOAT_T intensity, ///< the intensity of the peak that has been added -in
  FLOAT_T location ///< the location of the peak that has been added -in
  )
{
  // is new peak the smallest peak
  if(get_num_peaks() == 1 || 
     min_peak_mz > location){
    min_peak_mz = location;
  }
  // is new peak the largest peak
  if(get_num_peaks() == 1 || 
     max_peak_mz < location){
    max_peak_mz = location;
  }
  // update total_energy
  total_energy += intensity;
}


/**
 * \returns the number of the first scan
 */
int SPECTRUM_T::get_first_scan()
{
  return first_scan;
}

/**
 * sets the number of the first scan
 */
void SPECTRUM_T::set_first_scan(
  int first_scan ///< the first_scan -in
  )
{
  this->first_scan = first_scan;
}

/**
 * \returns the number of the last scan
 */
int SPECTRUM_T::get_last_scan(
  )
{
  return last_scan;
}

/**
 * sets the number of the last scan
 */
void SPECTRUM_T::set_last_scan(
  int last_scan ///< the last scan -in
  )
{
  this->last_scan = last_scan;
}

/**
 * \returns the spectrum_id
 */
int SPECTRUM_T::get_id(
  )
{
  return id;
}

/**
 * sets the spectrum_id
 */
void SPECTRUM_T::set_id(
  int id ///< the id -in
  )
{
  this->id = id;
}

/**
 * \returns the spectrum_type
 */
SPECTRUM_TYPE_T SPECTRUM_T::get_spectrum_type()
{
  return spectrum_type;
}

/**
 * sets the spectrum_type
 */
void SPECTRUM_T::set_spectrum_type(
  SPECTRUM_TYPE_T spectrum_type ///< the spectrum type -in
  )
{
  this->spectrum_type = spectrum_type;
}

/**
 * \returns the m/z of the precursor
 */
FLOAT_T SPECTRUM_T::get_precursor_mz()
{
  return precursor_mz;
}

/**
 * sets the m/z of the precursor
 */
void SPECTRUM_T::set_precursor_mz(
  FLOAT_T precursor_mz ///< the precursor_mz -in
  )
{
  this->precursor_mz = precursor_mz;
}

/**
 * \returns the minimum m/z of all peaks
 */
FLOAT_T SPECTRUM_T::get_min_peak_mz(
  )
{
  return min_peak_mz;
}

/**
 * \returns the maximum m/z of all peaks
 */
FLOAT_T SPECTRUM_T::get_max_peak_mz()
{
  return max_peak_mz;
}

/**
 * \returns the number of peaks
 */
int SPECTRUM_T::get_num_peaks(
  )
{
  return peaks.size();
}

/**
 * \returns the sum of intensities in all peaks
 */
double SPECTRUM_T::get_total_energy()
{
  return total_energy;
}


/**
 * \returns the filename of the ms2 file the spectrum was parsed
 * returns a char* to a heap allocated copy of the filename
 * user must free the memory
 */
char* SPECTRUM_T::get_filename()
{
  
  int filename_length = strlen(filename) +1; // +\0
  char * copy_filename = 
    (char *)mymalloc(sizeof(char)*filename_length);
  return strncpy(copy_filename, filename, filename_length);  
}

/**
 * sets the filename of the ms2 file the spectrum was parsed
 * copies the value from arguement char* filename into a heap allocated memory
 * frees memory for the filename that is replaced
 */
void SPECTRUM_T::set_filename(
  char* filename ///< the filename -in
  )
{
  free(filename);
  set_new_filename(filename);
}

/**
 * sets the filename of the ms2 file the spectrum was parsed
 * this function should be used only the first time the filename is set
 * to change existing filename use set_spectrum_filename
 * copies the value from arguement char* filename into a heap allocated memory
 */
void SPECTRUM_T::set_new_filename(
  char* filename ///< the filename -in
  )
{
  int filename_length = strlen(filename) +1; // +\0
  char * copy_filename = 
    (char *)mymalloc(sizeof(char)*filename_length);

  this -> filename =
    strncpy(copy_filename,filename,filename_length);  
}

/**
 * \returns the number of possible charge states of this spectrum
 */
int SPECTRUM_T::get_num_possible_z()
{
  return num_possible_z;
}

/**
 * \returns the possible charge states of this spectrum
 * returns an int* to a heap allocated copy of the src spectrum
 * thus, the user must free the memory
 * number of possible charge states can be gained by 
 * the get_num_possible_z function.
 */
int* SPECTRUM_T::get_possible_z()
{
  int num_possible_z_index = 0;
  int* new_possible_z = 
    (int*)mymalloc(sizeof(int)*num_possible_z);
  
  for(; num_possible_z_index < num_possible_z; 
      ++num_possible_z_index){
  
    new_possible_z[num_possible_z_index]
      = possible_z[num_possible_z_index];
  }
  return new_possible_z;
}

/**
 * \returns a pointer to an array of the possible charge states of this spectrum
 * User must NOT free this or alter, not a copy
 * number of possible charge states can be gained by 
 * the get_num_possible_z function.
 */
int* SPECTRUM_T::get_possible_z_pointer()
{
  return possible_z;
}

/**
 *  Considers all parameters and allocates an array of 
 *  charges that should be searched.  Returns the number of charges
 *  in that array.  Protects get_spectrum_possible_z_pointer, could make it
 *  private.
 */ 
int SPECTRUM_T::get_charges_to_search(int** select_charge_array){

  int total_charges = num_possible_z;
  int* all_charge_array = possible_z;

  int param_charge = 0;
  char* charge_str = get_string_parameter_pointer("spectrum-charge");
  int i=0;

  // Return full array of charges
  if( strcmp( charge_str, "all") == 0){

    *select_charge_array = (int*)mymalloc(sizeof(int) * total_charges);
    for(i=0; i < total_charges; i++){
      (*select_charge_array)[i] = all_charge_array[i];
    }
    return total_charges;
  }
  // else return one charge

  param_charge = atoi(charge_str);

  if( (param_charge < 1) || (param_charge > 3) ){
    carp(CARP_FATAL, "spectrum-charge option must be 1,2,3, or 'all'.  " \
         "%s is not valid", charge_str);
  }


  for(i=0; i<total_charges; i++){
     
    if( all_charge_array[i] == param_charge ){ 
      *select_charge_array = (int*)mymalloc(sizeof(int));
      **select_charge_array = param_charge;
      return 1;
    }
  }

  // Else none to be searched
  *select_charge_array = NULL;
  return 0;

}
/**
 * sets the possible charge states of this spectrum
 * this function should only be used when possible_z is set to NULL
 * to change existing possible_z use set_spectrum_possible_z()
 * the function copies the possible_z into a heap allocated memory
 * num_possible_z must match the array size of possible_z 
 * updates the number of possible charge states field
 */
void SPECTRUM_T::set_new_possible_z(
  int* possible_z, ///< possible z array -in
  int num_possible_z ///< possible z array size -in
  )
{
  
  int possible_z_index = 0;
  int* new_possible_z = 
    (int*)mymalloc(sizeof(int)*num_possible_z);
  
  for(; possible_z_index < num_possible_z; ++possible_z_index){
    new_possible_z[possible_z_index] = possible_z[possible_z_index];
  }
  
  possible_z = new_possible_z;
  num_possible_z = num_possible_z;

}

/**
 * sets the possible charge states of this spectrum
 * the function copies the possible_z into a heap allocated memory
 * num_possible_z must match the array size of possible_z 
 * frees the memory of the possible_z that is replaced
 * updates the number of possible charge states field
 */
void SPECTRUM_T::set_possible_z(
  int* possible_z, ///< possible z array -in
  int num_possible_z ///< possible z array size -in
  )
{
  free(possible_z);
  set_new_possible_z(possible_z, num_possible_z);
}

/**
 * \returns The intensity of the peak with the maximum intensity.
 */
FLOAT_T SPECTRUM_T::get_max_peak_intensity()
{
  int num_peak_index = 0;
  FLOAT_T max_intensity = -1;

  for(; num_peak_index < get_num_peaks(); ++num_peak_index){
    if(max_intensity <= peaks[num_peak_index].get_intensity()) {
      max_intensity = peaks[num_peak_index].get_intensity();
    }
  }
  return max_intensity; 
}


/**
 * \returns The mass of the charged precursor ion, according to the formula 
 * mass = m/z * charge
 */
FLOAT_T SPECTRUM_T::get_mass(
  int charge ///< the charge of precursor ion -in
  )
{
  return get_precursor_mz()*charge;
}

/**
 * \returns The mass of the neutral precursor ion, according to the formula 
 * mass = m/z * charge - mass_H * charge
 */
FLOAT_T SPECTRUM_T::get_neutral_mass(
  int charge ///< the charge of precursor ion -in
  )
{
  return (get_mass(charge) - MASS_H*charge); // TESTME
}

/**
 * \returns The mass of the singly charged precursor ion, according to the formula 
 * mass = m/z * charge - (mass_H * (charge - 1))
 */
FLOAT_T SPECTRUM_T::get_singly_charged_mass(
  int charge ///< the charge of the precursor ion -in
  )
{
  return (get_mass(charge) - MASS_H*(charge-1));  // TESTME
}


/**
 * serialize the spectrum in binary
 * Form,
 * <int: first_scan><int: last_scan><int: id><SPECTRUM_TYPE_T: spectrum_type>
 * <float: precursor_mz><float: retention_time>
 */
void SPECTRUM_T::serialize(
  FILE* file ///< output stream -out
  )
{
  // serialize the spectrum struct
  fwrite(this, sizeof(SPECTRUM_T), 1, file);
}

/**
 * Parse the spectrum from the serialized spectrum
 *\returns the parsed spectrum , else returns NULL for failed parse
 */
SPECTRUM_T* SPECTRUM_T::parse_spectrum_binary(
  FILE* file ///< output stream -out
  )
{
  carp(CARP_INFO,"creating a new spectrum object:%d",sizeof(SPECTRUM_T));
  SPECTRUM_T* spectrum = new SPECTRUM_T();
  carp(CARP_INFO,"done creating a new spectrum object:%d",sizeof(SPECTRUM_T));
  // get spectrum struct
  if(fread(spectrum, (sizeof(SPECTRUM_T)), 1, file) != 1){
    carp(CARP_ERROR, "serialized file corrupted, incorrect spectrum format");
    free(spectrum);
    return NULL;
  }
  
  spectrum->has_peaks = FALSE;
  
  return spectrum;
}

/***********************************************************************
 * Normalize peak intensities so that they sum to unity.
 ***********************************************************************/
void SPECTRUM_T::sum_normalize_spectrum()
{
  for (PEAK_ITERATOR_T peak_iterator = begin();
       peak_iterator != end();
       ++peak_iterator) {
    FLOAT_T new_intensity = peak_iterator -> get_intensity() / this->total_energy;
    peak_iterator -> set_intensity(new_intensity);
  }
}

/***********************************************************************
 * Populate peaks with rank information.
 ***********************************************************************/
void SPECTRUM_T::spectrum_rank_peaks()
{
  PEAK_ITERATOR_T peak_iterator;
  PEAK_T::sort_peaks(peaks, _PEAK_INTENSITY);
  sorted_by_intensity = TRUE;
  sorted_by_mz = FALSE;
  int rank = get_num_peaks();
  int num_peaks = get_num_peaks();

  for (peak_iterator = begin();
       peak_iterator != end();
       ++peak_iterator) {
    FLOAT_T new_rank = rank/(float)num_peaks;
    rank--;
    peak_iterator -> set_intensity_rank(new_rank);
  }
}



PEAK_ITERATOR_T SPECTRUM_T::begin() {
  return peaks.begin();
}

PEAK_ITERATOR_T SPECTRUM_T::end() {
  return peaks.end();
}




/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

