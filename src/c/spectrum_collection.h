/**
 * \file spectrum_collection.h 
 * AUTHOR: Chris Park
 * CREATE DATE: 28 June 2006
 * $Revision: 1.28 $
 * \brief Object for representing many spectra.
 *****************************************************************************/
#ifndef SPECTRUM_COLLECTION_H
#define SPECTRUM_COLLECTION_H

#include <stdio.h>
#include "objects.h"
#include "spectrum.h"
#include "carp.h"

#define MAX_SPECTRA 40000 ///< max number of spectrums
#define MAX_COMMENT 1000 ///< max length of comment

#ifdef __cplusplus
class SPECTRUM_COLLECTION_T {
 protected:

  int  num_spectra;     ///< The number of spectra
  int  num_charged_spectra; ///< The number of spectra assuming differnt charge(i.e. one spectrum with two charge states are counted as two spectra)
  char* filename;     ///< Optional filename
  BOOLEAN_T is_parsed; ///< Have we parsed all the spectra from the file?

 public:
  SPECTRUM_T* spectra[MAX_SPECTRA];  ///< The spectrum peaks
  char comment[MAX_COMMENT];    ///< The spectrum_collection header lines

  void init();
  SPECTRUM_COLLECTION_T();
  SPECTRUM_COLLECTION_T(char* filename);
  ~SPECTRUM_COLLECTION_T();
  
  void print(FILE* file);
  static void copy(SPECTRUM_COLLECTION_T* src,
		   SPECTRUM_COLLECTION_T* dest);
  BOOLEAN_T parse();
  
  BOOLEAN_T get_spectrum(int first_scan,  ///< The first scan of the spectrum to retrieve -in
			 SPECTRUM_T* spectrum
  );
 
  BOOLEAN_T add(SPECTRUM_T* spectrum);
  BOOLEAN_T add_to_end(SPECTRUM_T* spectrum);
  BOOLEAN_T remove(SPECTRUM_T* spectrum);
  
  void set_new_filename(char* filename);
  void set_filename(char* filename);
  
  char* get_filename();

  int get_num_spectra();

  int get_num_charged_spectra();
  
  char* get_comment();

  void set_comment(char* new_comment);

  BOOLEAN_T get_is_parsed();

  FILE** get_psm_result_filenames(
  char* psm_result_folder_name, ///< the folder name for where the result file should be placed -in
  char*** psm_result_filenames, ///< pointer to be set to the array of filenames for both the target and decoy psm results -out
  int number_decoy_set,  ///< the number of decoy sets to produce -in
  char* file_extension ///< the file extension of the spectrum file(i.e. ".ms2") -in
  );

  BOOLEAN_T serailize_header(
			     char* fasta_file,
			     FILE* psm_file);
  
  static BOOLEAN_T serialize_total_number_of_spectra(
					      int spectra_idx,
					      FILE* psm_file
					      );
  
  
};

#endif

#ifdef __cplusplus
extern "C" {
#endif


/******************************************************************************/

/**
 * Instantiates a new spectrum_iterator object from spectrum_collection.
 * \returns a SPECTRUM_ITERATOR_T object.
 */
SPECTRUM_ITERATOR_T* new_spectrum_iterator(
  SPECTRUM_COLLECTION_T* spectrum_collection///< spectrum_collection to iterate -in
);        

/**
 * Instantiates a new spectrum_iterator object from
 * spectrum_collection.  This iterator returns unique spectrum-charge
 * pairs (e.g.a spectrum to be searched as +2 and +3 is returned once
 * as +2 and once as +3).  The charge is returned by setting the int
 * pointer in the argument list.  The iterator also filters spectra by
 * mass so that none outside the spectrum-min-mass--spectrum-max-mass
 * range (as defined in parameter.c).
 * \returns a SPECTRUM_ITERATOR_T object.
 */
FILTERED_SPECTRUM_CHARGE_ITERATOR_T* new_filtered_spectrum_charge_iterator(
  SPECTRUM_COLLECTION_T* spectrum_collection///< spectra to iterate over
);        

/**
 * Frees an allocated spectrum_iterator object.
 */
void free_spectrum_iterator(
  SPECTRUM_ITERATOR_T* spectrum_iterator///< free spectrum_iterator -in
);

/**
 * The basic iterator function has_next.
 */
BOOLEAN_T spectrum_iterator_has_next(
  SPECTRUM_ITERATOR_T* spectrum_iterator///< is there a next spectrum? -in
);

/**
 * The basic iterator function has_next.
 */
BOOLEAN_T filtered_spectrum_charge_iterator_has_next(
  FILTERED_SPECTRUM_CHARGE_ITERATOR_T* iterator);

/**
 * The basic iterator function next.
 */
SPECTRUM_T* spectrum_iterator_next(
  SPECTRUM_ITERATOR_T* spectrum_iterator///< return the next spectrum -in
);

/**
 * The basic iterator function next.  Also returns the charge state to
 * use for this spectrum.
 */
SPECTRUM_T* filtered_spectrum_charge_iterator_next(
  FILTERED_SPECTRUM_CHARGE_ITERATOR_T* iterator,///< return spec from here -in
  int* charge);                 ///< put charge here -out


#ifdef __cplusplus
}
#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
