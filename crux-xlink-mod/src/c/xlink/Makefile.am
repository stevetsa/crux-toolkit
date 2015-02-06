# Automake template for xlink-search library
AM_CXXFLAGS = -Wall -W  $(USEDOUBLES) $(USEDARWIN) $(SEARCH_ENABLED) \
            -D_LARGEFILE_SOURCE=1 -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE


libxlink_la_SOURCES = \
  XHHC_Peptide.cpp \
  XHHC_Peptide.h \
  LinkedPeptide.cpp \
  LinkedPeptide.h \
  xhhc.cpp \
  xhhc.h \
  LinkedIonSeries.cpp \
  LinkedIonSeries.h \
  xhhc_score_peptide_spectrum.cpp \
  xhhc_scorer.cpp \
  xhhc_score_peptide_spectrum.cpp \
  xhhc_scorer.cpp \
  xhhc_scorer.h \
  xhhc_search.cpp \
  LinearPeptide.cpp \
  LinearPeptide.h \
  SelfLoopPeptide.cpp \
  SelfLoopPeptide.h \
  XLink.cpp \
  XLink.h \
  XLinkBondMap.cpp \
  XLinkBondMap.h \
  XLinkPeptide.cpp \
  XLinkPeptide.h \
  XLinkMatch.h \
  XLinkMatch.cpp \
  XLinkMatchCollection.h \
  XLinkMatchCollection.cpp \
  XLinkablePeptide.cpp \
  XLinkablePeptide.h \
  XLinkablePeptideIterator.cpp \
  XLinkablePeptideIterator.h \
  SearchForXLinks.cpp \
  SearchForXLinks.h \
  XLinkScorer.cpp \
  XLinkScorer.h \
  xlink_search.cpp \
  xlink_compute_qvalues.cpp \
  generate_xlink_peptides.cpp \
  XLinkSite.cpp \
  XLinkSite.h \
  XLinkBondMap.cpp \
  XLinkBondMap.h \
  XLinkablePeptide.cpp \
  XLinkablePeptide.h





libxlink_la_CXXFLAGS = $(AM_CXXFLAGS) \
  -I$(top_srcdir)/src/c/ \
  -I$(top_srcdir)/src/c/qranker-barista \
  -I$(top_srcdir)/src/external/MSToolkit/include

noinst_LTLIBRARIES = libxlink.la

