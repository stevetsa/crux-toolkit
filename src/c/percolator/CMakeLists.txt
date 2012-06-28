project(percolator)

cmake_minimum_required(VERSION 2.8.1)
cmake_policy(VERSION 2.8.1)

set_property(
  DIRECTORY
  PROPERTY
  COMPILE_DEFINITIONS
  _LARGEFILE_SOURCE
)
include_directories(..)
include_directories(../ArrayLibrary)
add_library(
  percolator
  STATIC
  BaseSpline.cpp
  Caller.cpp
  DataSet.cpp
  Globals.cpp
  IntraSetRelation.cpp
  LogisticRegression.cpp
  Normalizer.cpp
  Option.cpp
  PercolatorCInterface.cpp
  PosteriorEstimator.cpp
  ResultHolder.cpp
  SanityCheck.cpp
  Scores.cpp
  SetHandler.cpp
  SqtSanityCheck.cpp
  StdvNormalizer.cpp
  UniNormalizer.cpp
  ssl.cpp
)
