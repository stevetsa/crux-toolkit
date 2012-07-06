cmake_minimum_required(VERSION 2.8.1)
cmake_policy(VERSION 2.8.1)

add_subdirectory(ArrayLibrary)
add_subdirectory(bullseye)
add_subdirectory(hardklor)
add_subdirectory(percolator)
add_subdirectory(posterior-error)
add_subdirectory(qranker-barista)
add_subdirectory(xlink)

include_directories(.)
include_directories(ArrayLibrary)
include_directories(bullseye)
include_directories(hardklor)
include_directories(percolator)
include_directories(qranker-barista)
include_directories(posterior-error)
include_directories(xlink)
include_directories(../external)
include_directories(../external/MSToolkit/include)
link_directories(../external/MSToolkit)
link_directories(.)
link_directories(./ArrayLibrary)
link_directories(./bullseye)
link_directories(./hardklor)
link_directories(./percolator)
link_directories(./posterior-error)
link_directories(./qranker-barista)
link_directories(./xlink)

set (
  crux_lib_files
  Alphabet.cpp
  analyze_psms.cpp
  carp.cpp
  ComputeQValues.cpp
  CreateIndex.cpp
  crux-utils.cpp
  Database.cpp
  DatabaseProteinIterator.cpp
  DatabasePeptideIterator.cpp
  CruxApplication.cpp
  CruxApplicationList.cpp
  DelimitedFile.cpp
  DelimitedFileReader.cpp
  DelimitedFileWriter.cpp
  ExtractColumns.cpp
  ExtractRows.cpp
  FilteredSpectrumChargeIterator.cpp
  GeneratePeptidesIterator.cpp
  GetMs2Spectrum.cpp
  hash.cpp
  Index.cpp
  IndexFile.cpp
  IndexPeptideIterator.cpp
  IndexFilteredPeptideIterator.cpp
  Ion.cpp
  IonConstraint.cpp
  IonFilteredIterator.cpp
  IonSeries.cpp
  linked_list.cpp
  mass.cpp
  Match.cpp
  MatchColumns.cpp
  MatchFileReader.cpp
  MatchFileWriter.cpp
  MatchSearch.cpp
  MatchCollection.cpp
  MatchCollectionIterator.cpp
  MatchIterator.cpp
  MGFSpectrumCollection.cpp
  modifications.cpp
  ModifiedPeptidesIterator.cpp
  MS2SpectrumCollection.cpp
  MSToolkitSpectrumCollection.cpp
  OutputFiles.cpp
  parameter.cpp
  parse_arguments.cpp
  Peak.cpp
  Peptide.cpp
  PeptideConstraint.cpp
  peptide_modifications.cpp
  PeptideSrc.cpp
  PepXMLWriter.cpp
  Percolator.cpp
  GeneratePeptides.cpp
  PredictPeptideIons.cpp
  PrintProcessedSpectra.cpp
  PrintVersion.cpp
  Protein.cpp
  ProteinPeptideIterator.cpp
  ProteinIndex.cpp
  ProteinIndexIterator.cpp
  q-value.cpp
  SequestSearch.cpp
  SortColumn.cpp
  Scorer.cpp
  sorter.cpp
  SpectralCounts.cpp
  SpectrumCollection.cpp
  SpectrumCollectionFactory.cpp
  Spectrum.cpp
  SpectrumZState.cpp
  StatColumn.cpp
  utils.cpp
)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set(
    crux_lib_files
    ${crux_lib_files}
    WinCrux.cpp
  )
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

add_library(
  crux-support
  STATIC
  ${crux_lib_files}
)

set_property(
  SOURCE
  crux-main.cpp
  PROPERTY
  COMPILE_DEFINITIONS
  MAIN
)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_property(
      DIRECTORY
      PROPERTY
      COMPILE_DEFINITIONS
      _CRT_SECURE_NO_WARNINGS
      SEARCH_ENABLED
    )
else (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_property(
	  DIRECTORY
	  PROPERTY
	  COMPILE_DEFINITIONS
	  _LARGEFILE_SOURCE
	  _FILE_OFFSET_BITS=64
	  _GNU_SOURCE
	  GCC
	  HAVE_EXPAT_CONFIG_H
    SEARCH_ENABLED
	)
	if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	  set_property(
		DIRECTORY
		APPEND
		PROPERTY
		COMPILE_DEFINITIONS
		DARWIN
	  )
	  if (${CMAKE_SYSTEM_VERSION} VERSION_LESS "11")
		set_property(
		  DIRECTORY
		  PROPERTY
		  COMPILE_DEFINITIONS
		  PRE_LION
		)
       endif(${CMAKE_SYSTEM_VERSION} VERSION_LESS "11")
    endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

add_executable(
  crux
  crux-main.cpp
)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
target_link_libraries(
  crux
  barista
  bullseye
  crux-support
  hardklor
  MSToolkit
  percolator
  posterior_error
  xlink
  ws2_32
  ArrayLibrary
)
else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
target_link_libraries(
  crux
  barista
  bullseye
  crux-support
  dl
  hardklor
  m
  MSToolkit
  percolator
  posterior_error
  pthread
  xlink
  ArrayLibrary
)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
