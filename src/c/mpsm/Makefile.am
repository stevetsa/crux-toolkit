# Automake template for mpsm-search library
libmpsm_la_SOURCES = \
  ChargeIndex.cpp \
  ChargeIndex.h \
  MPSM_ChargeMap.cpp \
  MPSM_ChargeMap.h \
  MPSM_Match.cpp \
  MPSM_Match.h \
  MPSM_MatchCollection.cpp \
  MPSM_MatchCollection.h \
  MPSM_OutputFiles.cpp \
  MPSM_OutputFiles.h \
  MPSM_Scorer.cpp \
  MPSM_Scorer.h \
  RetentionPredictor.cpp \
  RetentionPredictor.h \
  AKlammerRetentionPredictor.cpp \
  AKlammerRetentionPredictor.h \
  KrokhinRetentionPredictor.h \
  KrokhinRetentionPredictor.cpp \
  PalmbaldRetentionPredictor.cpp \
  PalmbaldRetentionPredictor.h \
  BigSmallReal.cpp \
  BigSmallReal.h \
  mpsm_search.cpp

libmpsm_la_CXXFLAGS = $(AM_CXXFLAGS) \
  -D SEARCH_ENABLED \
  -I$(top_srcdir)/src/c/ \
  -I$(top_srcdir)/src/c/MSToolkit \
  -I$(top_srcdir)/src/c/MSToolkit/zLib \
  -I$(top_srcdir)/src/c/MSToolkit/RAMP

noinst_LTLIBRARIES = libmpsm.la

