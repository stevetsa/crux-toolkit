/*************************************************************************//**
 * \file MatchColumns.h
 * \brief Just keeps track of column names for match files.
 **************************************************************************/

#ifndef MATCHCOLUMNS_H
#define MATCHCOLUMNS_H

//#define NEW_COLUMNS 1

enum MATCH_COLUMNS_T {
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
#ifdef NEW_COLUMNS
  WEIBULL_QVALUE_COL,
  WEIBULL_PEPTIDE_QVALUE_COL,      // NEW
  DECOY_XCORR_QVALUE_COL,
  DECOY_XCORR_PEPTIDE_QVALUE_COL,  // NEW
  PERCOLATOR_SCORE_COL,
  PERCOLATOR_RANK_COL,
  PERCOLATOR_QVALUE_COL,
  PERCOLATOR_PEPTIDE_QVALUE_COL,   // NEW
  QRANKER_SCORE_COL,
  QRANKER_QVALUE_COL,
  QRANKER_PEPTIDE_QVALUE_COL,      // NEW
#else
  WEIBULL_QVALUE_COL,
  DECOY_XCORR_QVALUE_COL,
  PERCOLATOR_SCORE_COL,
  PERCOLATOR_RANK_COL,
  PERCOLATOR_QVALUE_COL,
  QRANKER_SCORE_COL,
  QRANKER_QVALUE_COL,
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
  NZSTATE_COL,
  RTIME_MAX_DIFF_COL,
  PEPTIDES_SPECTRUM_COL,
  XCORR_SUM_DIFF_COL,
  XCORR_MAX_DIFF_COL,
  MZ1_AREA_RATIO_COL,
  SIN_SCORE_COL,
  NSAF_SCORE_COL,
  PARSIMONY_RANK_COL,

  NUMBER_MATCH_COLUMNS,
  INVALID_COL
};

/**
 * Get the name of a given column, by index.
 */
const char* get_column_header(
  int columnIndex
);

#endif // MATCHCOLUMNS_H
