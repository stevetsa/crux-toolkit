project(barista)

cmake_minimum_required(VERSION 2.8.1)
cmake_policy(VERSION 2.8.1)

include_directories(..)

set_property(
  DIRECTORY
  PROPERTY
  COMPILE_DEFINITIONS
  CRUX
  SEARCH_ENABLED
)
if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_property(
	  DIRECTORY
	  PROPERTY
	  COMPILE_DEFINITIONS
      _CRT_SECURE_NO_WARNINGS
      _SCL_SECURE_NO_WARNINGS
    )
else (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_property(
	  DIRECTORY
	  PROPERTY
	  COMPILE_DEFINITIONS
	  _FILE_OFFSET_BITS=64
	  _GNU_SOURCE
	  _LARGEFILE_SOURCE
	  GCC
	)
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

add_library(
  barista
  STATIC
  Barista.cpp
  ProtScores.cpp
  NeuralNet.cpp
  SpecFeatures.cpp
  BipartiteGraph.cpp
  PSMScores.cpp
  PepScores.cpp
  SQTParser.cpp
  TabDelimParser.cpp
  DataSet.cpp
  QRanker.cpp
)
