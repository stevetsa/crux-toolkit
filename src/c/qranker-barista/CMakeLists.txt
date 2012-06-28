project(barista)

cmake_minimum_required(VERSION 2.8.1)
cmake_policy(VERSION 2.8.1)

set_property(
  DIRECTORY
  PROPERTY
  COMPILE_DEFINITIONS
  _LARGEFILE_SOURCE
  _FILE_OFFSET_BITS=64
  _GNU_SOURCE
)
include_directories(..)
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
