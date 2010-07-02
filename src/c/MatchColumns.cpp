/*************************************************************************//**
 * \file MatchColumns.cpp
 * \brief Just keeps track of column names for match files.
 ****************************************************************************/

#include "MatchColumns.h"

static const char* match_column_strings[NUMBER_MATCH_COLUMNS] = {
  "scan",
  "charge",
  "spectrum precursor m/z",
  "spectrum neutral mass",
  "peptide mass",
  "delta_cn",
  "sp score",
  "sp rank",
  "xcorr score",
  "xcorr rank",
  "p-value",
 #ifdef NEW_COLUMNS
  "Weibull PSM q-value",
  "Weibull peptide q-value",
  "decoy PSM q-value",
  "decoy peptide q-value",
  "percolator score",
  "percolator rank",
  "percolator PSM q-value",
  "percolator peptide q-value",
  "q-ranker score",
  "q-ranker PSM q-value",
  "q-ranker peptide q-value",
#else
  "Weibull est. q-value",
  "decoy q-value (xcorr)",
  "percolator score",
  "percolator rank",
  "percolator q-value",
  "q-ranker score",
  "q-ranker q-value",
#endif
  "b/y ions matched",
  "b/y ions total",
  "matches/spectrum",
  "sequence",
  "cleavage type",
  "protein id",
  "flanking aa",
  "unshuffled sequence",
  "eta",
  "beta",
  "shift",
  "corr"
};

/**
 * Get the name of a given column, by index.
 */
const char* get_column_header(
  int columnIndex
) {
  return(match_column_strings[columnIndex]);
}
