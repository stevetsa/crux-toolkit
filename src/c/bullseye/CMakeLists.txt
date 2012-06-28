project(bullseye)

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
include_directories(../hardklor)
include_directories(../../external/MSToolkit/include)
add_library(
  bullseye
  STATIC
	bullseye.cpp
	CKronik2.cpp
  CruxBullseyeApplication.cpp
)
