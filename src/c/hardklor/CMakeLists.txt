project(hardklor)

cmake_minimum_required(VERSION 2.8.1)
cmake_policy(VERSION 2.8.1)

set_property(
  DIRECTORY
  PROPERTY
  COMPILE_DEFINITIONS
  GCC
  _LARGEFILE_SOURCE
  HAVE_EXPAT_CONFIG_H
  _FILE_OFFSET_BITS=64
  _GNU_SOURCE
  CRUX
)
include_directories(..)
include_directories(../../external/MSToolkit/include)
add_library(
  hardklor
  STATIC
  CAveragine.cpp
  CHardklor.cpp
  CHardklor2.cpp
  CHardklorFileReader.cpp
  CHardklorParser.cpp
  CHardklorProtein.cpp
  CHardklorSetting.cpp
  CHardklorVariant.cpp
  CMercury8.cpp
  CModelLibrary.cpp
  CNoiseReduction.cpp
  CPeriodicTable.cpp
  CSettings.cpp
  CSpecAnalyze.cpp
  CSplitSpectrum.cpp
  FFT.cpp
  FFT-HK.cpp
  HardklorApp.cpp
  S2N.cpp
  Smooth.cpp
  SpecAnalyzeSupport.cpp
  CruxHardklorApplication.cpp
)
