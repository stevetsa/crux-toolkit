/*************************************************************************//**
 * \file MatchColumns.h
 * \brief Just keeps track of column names for match files.
 **************************************************************************/

#ifndef MATCHCOLUMNS_H
#define MATCHCOLUMNS_H

enum _match_columns {
  SCAN_COL,
  CHARGE_COL,
  SPECTRUM_PRECURSOR_MZ_COL,
  SPECTRUM_NEUTRAL_MASS_COL,
  PEPTIDE_MASS_COL,
  DELTA_CN_COL,
  SP_SCORE_COL,
  SP_RANK_COL,
  XCORR_SCORE_COL,
  XCORR_RANK_COL,
  PVALUE_COL,
  WEIBULL_QVALUE_COL,
#ifdef NEW_COLUMNS
  WEIBULL_PEPTIDE_QVALUE_COL,
#endif
  DECOY_XCORR_QVALUE_COL,
#ifdef NEW_COLUMNS
  DECOY_XCORR_PEPTIDE_QVALUE_COL,
#endif
  PERCOLATOR_SCORE_COL,
  PERCOLATOR_RANK_COL,
  PERCOLATOR_QVALUE_COL,
#ifdef NEW_COLUMNS
  PERCOLATOR_PEPTIDE_QVALUE_COL,
#endif
  QRANKER_SCORE_COL,
  QRANKER_QVALUE_COL,
#ifdef NEW_COLUMNS
  QRANKER_PEPTIDE_QVALUE_COL,
#endif
  BY_IONS_MATCHED_COL,
  BY_IONS_TOTAL_COL,
  MATCHES_SPECTRUM_COL,
  SEQUENCE_COL,
  CLEAVAGE_TYPE_COL,
  PROTEIN_ID_COL,
  FLANKING_AA_COL,
  UNSHUFFLED_SEQUENCE_COL,
  ETA_COL,
  BETA_COL,
  SHIFT_COL,
  CORR_COL,
  NUMBER_MATCH_COLUMNS
};

typedef enum _match_columns MATCH_COLUMNS_T;

/**
 * Get the name of a given column, by index.
 */
const char* get_column_header(
  int columnIndex
);

#endif // MATCHCOLUMNS_H
