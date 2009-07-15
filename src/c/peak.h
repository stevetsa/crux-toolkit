/**
 * \file peak.h
 * $Revision: 1.14 $
 * \brief Object for representing one peak in a spectrum.
 *
 * A peak is primarily identified via its intensity (height) and location
 * (position on the m/z axis).
 *
 */
#ifndef PEAK_H
#define PEAK_H

#ifdef __cplusplus
#include <vector>
#endif

#include <stdio.h>
#include <stdlib.h>
#include "objects.h"
#include "utils.h"

#ifdef __cplusplus
class PEAK_T{
 protected:
  FLOAT_T intensity;  ///< The intensity of the peak.
  FLOAT_T intensity_rank;  ///< The rank intensity of the peak.
  FLOAT_T location;   ///< The location of the peak.
 public:

  void init();
  PEAK_T();
  

  /**
   * \returns A PEAK_T object
   */

  PEAK_T(
	 FLOAT_T intensity, ///< intensity for the new peak -in 
	 FLOAT_T location ///< location for the new peak -in
	 );
  /**
   * \frees A PEAK_T object
   */
  ~PEAK_T();


  /**
   * \returns the intensity of PEAK_T object
   */
  FLOAT_T get_intensity();

  /**
   * sets the intensity rank of PEAK_T object
   */
  FLOAT_T get_intensity_rank();

  /**
   * \returns the location of PEAK_T object
   */
  FLOAT_T get_location();

  /**
   * sets the intensity of PEAK_T object
   */
  void set_intensity(
		     FLOAT_T intensity ///< the intensity -in
		     );
  /**
   * sets the intensity rank of PEAK_T object
   */
  void set_intensity_rank(
			  FLOAT_T intensity_rank ///< the intensity -in
			  );

  /**
   * sets the location of PEAK_T object
   */
  void set_location(
			 FLOAT_T location ///< the location -in
			 );

  /**
   * \prints the intensity and location of PEAK_T object to stdout
   */
  void print_peak();

  /**
   * \returns A heap allocated PEAK_T object array
   */
  static PEAK_T* allocate_peak_array(
  int num_peaks///< number of peaks to allocate -in
  );

 
  /**
 * \frees A PEAK_T object array
 */
  static void free_peak_array(
				 PEAK_T* garbage_peak ///<the peak array to free -in
				 ); 

  /**
   *\returns a pointer to the peak in the peak_array
   */
  static PEAK_T* find_peak(
		    PEAK_T* peak_array,///< peak_array to search -in
		    int index ///< the index of peak to fine -in
		    );
  static PEAK_T* find_peak(
		  std::vector<PEAK_T>& peak_vector,
		  unsigned int index);

  
  /**
   * sort peaks by their intensity or location
   * use the lib. function, qsort()
   */
  static void sort_peaks(
			 PEAK_T* peak_array, ///< peak array to sort -in/out
			 int num_peaks,  ///< number of total peaks -in
			 PEAK_SORT_TYPE_T sort_type ///< the sort type(location or intensity)
			 );

  static void sort_peaks(
			 std::vector<PEAK_T>& peak_vector,
			 PEAK_SORT_TYPE_T sort_type);

  friend bool compare_peaks_by_intensity(const PEAK_T& peak_1,
				  const PEAK_T& peak_2);

  
  friend bool compare_peaks_by_mz(const PEAK_T& peak_1,
			   const PEAK_T& peak_2);
};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

//#ifdef __cplusplus
//}
//#endif


#endif
