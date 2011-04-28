# Automake template for mpsm-search library
libmpsm_la_SOURCES = \
  ZStateIndex.cpp \
  ZStateIndex.h \
  MPSM_Match.cpp \
  MPSM_Match.h \
  MPSM_MatchCollection.cpp \
  MPSM_MatchCollection.h \
  MPSM_OutputFiles.cpp \
  MPSM_OutputFiles.h \
  MPSM_Scorer.cpp \
  MPSM_Scorer.h \
  MPSM_ZStateMap.cpp \
  MPSM_ZStateMap.h \
  RetentionPredictor.cpp \
  RetentionPredictor.h \
  AKlammerRetentionPredictor.cpp \
  AKlammerRetentionPredictor.h \
  AKlammerStaticRetentionPredictor.cpp \
  AKlammerStaticRetentionPredictor.h \
  KrokhinRetentionPredictor.h \
  KrokhinRetentionPredictor.cpp \
  NullRetentionPredictor.h \
  NullRetentionPredictor.cpp \
  PalmbaldRetentionPredictor.cpp \
  PalmbaldRetentionPredictor.h \
  AssignIonsMPSMSpectrum.h \
  AssignIonsMPSMSpectrum.cpp \
  ScoreMPSMSpectrum.h \
  ScoreMPSMSpectrum.cpp \
  SearchForMPSMS.h \
  SearchForMPSMS.cpp \
  PredictMPSMIons.h \
  PredictMPSMIons.cpp\
  svm.h \
  svm.cpp



libmpsm_la_CXXFLAGS = $(AM_CXXFLAGS) \
  -D SEARCH_ENABLED \
  -I$(top_srcdir)/src/c/ \
  -I$(top_srcdir)/src/c/MSToolkit \
  -I$(top_srcdir)/src/c/MSToolkit/zLib \
  -I$(top_srcdir)/src/c/MSToolkit/RAMP

noinst_LTLIBRARIES = libmpsm.la
