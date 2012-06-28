project(posterior-error)

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
  posterior_error
  STATIC
  BaseSpline.cpp
  LogisticRegression.cpp
  Numerical.cpp
  PackedMatrix.cpp
  PackedVector.cpp
  PosteriorEstimator.cpp
  Set.cpp
  Vector.cpp
  Vector.h 
)
