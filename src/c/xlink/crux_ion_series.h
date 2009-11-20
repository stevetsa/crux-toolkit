#ifndef CRUX_ION_SERIES_H_
#define CRUX_ION_SERIES_H_

#ifdef _cplusplus
extern "C" {
#endif

#include "ion_series.h"

  /*********************************/
  /*EXPOSE INTERNAL CRUX FUNCTIONS */
  /*********************************/
/**
 * \brief Find instances of amino acid which can incur neutral
 * losses: H2O (S|T|E|D), NH3(R|K|Q|N).  
 * Set the count of those observed so far for each cleavage index.
 * If no instance of amino acid, the count is assigned to 0
 * The information is used to determine if how many nh3 or h2o neutral
 * losses are possible. 
 */
void scan_for_aa_for_neutral_loss(
  ION_SERIES_T* ion_series ///< ion_series to print -in/out
  );

/**
 * Creates all the ions with no modifications up to the max charge
 * Adds each ion to ion_series
 *\returns TRUE if successfully generates all the ions, else FALSE
 */
BOOLEAN_T generate_ions_no_modification(
  ION_SERIES_T* ion_series, ///< ion_series to modify -in/out
  FLOAT_T* mass_matrix ///< the mass matrix that stores the mass
  );

/**
 * creates all the ions with specific modifications up to the max charge
 * copies all the existing ions that can be modified,
 * then applies the different modifications then adds the new modified ions to ion_series
 *\returns TRUE if successfully generates all the ions with modifications, else FALSE
 */
BOOLEAN_T generate_ions(
  ION_SERIES_T* ion_series, ///< ion_series to print -in/out  
  ION_MODIFICATION_T mod_type ///< generate ions of this modification_type -in/out
  );

/**
 * creates all the flanking ions up to the max charge
 * can only create flanking ions that are B|Y ions and don't have modification
 * assumes the ions with no modification all are at the begining of the ion[] in ion_series
 * copies all the existing ions that can be modified,
 * then applies the different modifications then adds the new modified ions to ion_series
 *\returns TRUE if successfully generates all the ions with modifications, else FALSE
 */
BOOLEAN_T generate_ions_flank(
  ION_SERIES_T* ion_series ///< ion_series to print -in/out
  );
#ifdef _cplusplus
}
#endif
#endif
