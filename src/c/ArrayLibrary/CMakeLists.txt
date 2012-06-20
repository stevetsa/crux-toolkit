project(crux)

cmake_minimum_required(VERSION 2.6)

add_library(
  array
  STATIC
  MatrixFunctions.cpp
  Numerical.cpp
  PackedMatrix.cpp
  PackedSetFunctions.cpp
  PackedVec.cpp
  Random.cpp
  VectorFunctions.cpp
)
