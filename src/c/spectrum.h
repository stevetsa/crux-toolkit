/**
 * \file spectrum.h 
 * $Revision: 1.43 $
 * \brief Object for representing one spectrum.
 *****************************************************************************/
#ifndef SPECTRUM_H
#define SPECTRUM_H

#ifdef __cplusplus
//MSToolkit Includes
#include "Spectrum.h"
#endif

//#include <stdio.h>
#include "utils.h"
#include "objects.h"
#include "peak.h"



#define MAX_I_LINES 2 // number of 'I' lines albe to parse for one spectrum object
#define MAX_D_LINES 2 // number of 'D' lines albe to parse for one spectrum object

#ifdef __cplusplus

typedef std::vector<PEAK_T>::iterator PEAK_ITERATOR_T;


/**
 * \class SPECTRUM_T 
 * \brief A mass spectrum

 * A mass spectrum consists mainly of a list of peak objects along with
 * some identifying information. A single spectrum is generated from one 
 * or more "scans" of the mass spectrometer; each scan is identified by 
 * a unique increasing positive integer. The range of scans that
 * generated a particular spectrum are indicated by the member variables 
 * "first_scan" and "last_scan". In addition to scan information, 
 * a tandem fragmentation mass spectrum has information 
 * about the m/z of the intact ion that generated the spectrum, which is
 * indicated by "precursor_mz" member variable.
 * Also, while the m/z of particular spectrum is known, the charge state of
 * the originating ion is unknown; the possible charge states of the 
 * precursor ion is stored "possible_z" and "num_possible_z". 
 * Finally, some summary information that can be derived from the spectrum
 * peaks but is convenient to have is stored as "min_peak_mz",
 * "max_peak_mz", and "total_energy".
 */
class SPECTRUM_T {
 protected:
  //MSToolkit::Spectrum mst_spectrum;

  int              first_scan;    ///< The number of the first scan
  int              last_scan;     ///< The number of the last scan
  int              id;            ///< A unique identifier
                                  // FIXME, this field is not set when parsing
  SPECTRUM_TYPE_T  spectrum_type; ///< The type of spectrum. 
  FLOAT_T            precursor_mz;  ///< The m/z of precursor (MS-MS spectra)
  int*             possible_z;    ///< The possible charge states of this spectrum
  int              num_possible_z;///< The number of possible charge states of this spectrum

  FLOAT_T            min_peak_mz;   ///< The minimum m/z of all peaks
  FLOAT_T            max_peak_mz;   ///< The maximum m/z of all peaks
    //int              num_peaks;     ///< The number of peaks
  double           total_energy;  ///< The sum of intensities in all peaks
  char*            filename;      ///< Optional filename
  char*            i_lines[MAX_I_LINES]; ///< store i lines, upto MAX_I_LINES
  char*            d_lines[MAX_D_LINES]; ///< store d lines, upto MAX_D_LINES 
  BOOLEAN_T        sorted_by_mz; ///< Are the spectrum peaks sorted by m/z...
  BOOLEAN_T        sorted_by_intensity; ///< ... or by intensity?
  BOOLEAN_T        has_mz_peak_array; ///< Is the mz_peak_array populated.

  /*
   * Parses the 'S' line of a spectrum
   * \returns TRUE if success. FALSE is failure.
   */
  BOOLEAN_T parse_S_line(
    char* line, ///< 'S' line to parse -in
    int buf_length ///< line length -in
  );

  /*
   * Parses the 'Z' line of the a spectrum
   * \returns TRUE if success. FALSE is failure.
   */
  BOOLEAN_T parse_Z_line(
    char* line  ///< 'Z' line to parse -in
  );

/*
 * Parses the 'D' line of the a spectrum
 * \returns TRUE if success. FALSE is failure.
 */
BOOLEAN_T parse_D_line(
  char* line  ///< 'D' line to parse -in
);

/*
 * Parses the 'I' line of the a spectrum
 * \returns TRUE if success. FALSE is failure.
 */
BOOLEAN_T parse_I_line(
  char* line  ///< 'I' line to parse -in
);

/*
 * Adds a possible charge(z) to the spectrum.
 * Must not exceed the MAX_CHARGE capacity
 */
BOOLEAN_T add_possible_z(
  int charge  ///< charge to add
);

 void populate_mz_peak_array();



 public:
 //TODO : Protect these!
  BOOLEAN_T        has_peaks;  ///< Does the spectrum contain peak information
  PEAK_T**         mz_peak_array;  ///< Allows rapid peak retrieval by mz.
    
  vector<PEAK_T>   peaks;  
  //PEAK_T*          peaks;         ///< The spectrum peaks


  void init();
  SPECTRUM_T();
  SPECTRUM_T(
    int               first_scan,         ///< The number of the first scan -in
    int               last_scan,          ///< The number of the last scan -in
    SPECTRUM_TYPE_T   spectrum_type,      ///< The type of spectrum. -in
    FLOAT_T             precursor_mz,       ///< The m/z of the precursor (for MS-MS spectra) -in
    int*              possible_z,         ///< The possible charge states of this spectrum  -in
    int               num_possible_z,     ///< The number of possible charge states of this spectrum  -in  
    char*             filename);          ///< Optional filename  -in    
    

 SPECTRUM_T(MSToolkit::Spectrum& mst_spectrum);

 ~SPECTRUM_T(); // free_spectrum
 void print(FILE* file);
 void print_sqt(FILE* file,///< output file to print at -out
  int num_matches,      ///< number of peptides compared to this spec -in
  int charge            ///< charge used for the search -in
  ); 
  
 void print_stdout();

 static void copy(
   SPECTRUM_T* src, ///< the source spectrum -in
   SPECTRUM_T* dest ///< the destination spectrum -out
   );

 BOOLEAN_T parse_spectrum_file(
  FILE* file, ///< the input file stream -in
  char* filename ///< filename of the spectrum, should not free -in
  );

 BOOLEAN_T parse_spectrum(
  char*      filename ///< the file to parse -in
  );

 static SPECTRUM_T* parse_spectrum_binary(
   FILE* file ///< output stream -out
  );
  
 void sum_normalize_spectrum();
 void spectrum_rank_peaks();
 int get_first_scan();
 void set_first_scan(int first_scan);
 int get_last_scan();
 void set_last_scan(int last_scan);

 int get_id();
 void set_id(int id);
 
 SPECTRUM_TYPE_T get_spectrum_type();
 
 void set_spectrum_type(SPECTRUM_TYPE_T);


 FLOAT_T get_precursor_mz();
 void set_precursor_mz(FLOAT_T precursor_mz);
 
 

 int* get_possible_z();
 
 int* get_possible_z_pointer();
 
 int get_charged_to_search(int**);
 void set_possible_z(int* possible_z, int num_possible_z);
 
 void set_new_possible_z(
   int* possible_z, ///< possible z array -in
   int num_possible_z ///< possible z array size -in
 );		 

 int get_num_possible_z();

 FLOAT_T get_min_peak_mz();

 FLOAT_T get_max_peak_mz();
 
 int get_num_peaks();

 double get_total_energy();

 char* get_filename();

 void set_filename(char* filename);

 void set_new_filename(char* filename);
 
 FLOAT_T get_max_peak_intensity();

 FLOAT_T get_mass(int charge);
 
 FLOAT_T get_neutral_mass(int charge);
 
 FLOAT_T get_singly_charged_mass(int charge);

 void update_spectrum_fields(
   FLOAT_T intensity, ///< the intensity of the peak that has been added -in
   FLOAT_T location ///< the location of the peak that has been added -in
   );

 BOOLEAN_T add_peak(
  FLOAT_T intensity, ///< the intensity of peak to add -in
  FLOAT_T location_mz ///< the location of peak to add -in
  );
 PEAK_T* get_nearest_peak(
  FLOAT_T mz, ///< the mz of the peak around which to sum intensities -in
  FLOAT_T max ///< the maximum distance to get intensity -in
  );
 
 FLOAT_T get_nearby_intensity_sum( 
  FLOAT_T mz,             ///< the mz of the peak around which to sum intensities
  FLOAT_T tol             ///< the tolerance within which to sum intensities
  );

 
/**
 * process the spectrum, according the score type
 *\returns a new spectrum that has been preprocessed
 */
SPECTRUM_T* process_spectrum(
  SPECTRUM_T* spectrum, ///< the spectrum to processes -in
  SCORER_TYPE_T score_type ///< the score type to which the spectrum should be sorted -in
  );

/**
 * serialize the spectrum in binary
 * Form,
 * <int: first_scan><int: last_scan><int: id><SPECTRUM_TYPE_T: spectrum_type>
 * <float: precursor_mz><float: retention_time>
 */
void serialize(
  FILE* file ///< output stream -out
  );

/***********************************************************************
 * Normalize peak intensities so that they sum to unity.
 ***********************************************************************/
void sum_normalize_spectrum(
  SPECTRUM_T* spectrum
  );

/***********************************************************************
 * Populate peaks with rank information.
 ***********************************************************************/
void rank_peaks(
  SPECTRUM_T* spectrum
  );

 int get_charges_to_search(int** select_charge_array);  



 PEAK_ITERATOR_T begin();
 PEAK_ITERATOR_T end();

};
#endif



/**
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

#endif

/** \mainpage The crux API documentation page.
 * \section Introduction
 * Welcome to crux, a C software package for analysis of tandem mass
 * spectrometry data. Click on the links above to see documentation for
 * crux objects and their user interfaces.
 */
