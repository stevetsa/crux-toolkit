project(crux)

cmake_minimum_required(VERSION 2.6)

set_property(
  DIRECTORY
  PROPERTY
  COMPILE_DEFINITIONS
  SEARCH_ENABLED
  USEDARWIN
  _LARGEFILE_SOURCE
  _FILE_OFFSET_BITS=64
  _GNU_SOURCE
)
include_directories(..)
include_directories(../external/MSToolkit/include)
add_library(
  xlink
  STATIC
  XHHC_Peptide.cpp
  LinkedPeptide.cpp
  xhhc.cpp
  LinkedIonSeries.cpp
  xhhc_score_peptide_spectrum.cpp
  xhhc_scorer.cpp
  xhhc_score_peptide_spectrum.cpp
  xhhc_scorer.cpp
  xhhc_search.cpp
  xlink_search.cpp
  SearchForXLinks.cpp
  xlink_compute_qvalues.cpp
  XLinkSite.cpp
  XLinkBondMap.cpp
  XLinkablePeptide.cpp
)
