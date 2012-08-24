project(comet)

cmake_minimum_required(VERSION 2.8.1)
cmake_policy(VERSION 2.8.1)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_property(
      DIRECTORY
      PROPERTY
      COMPILE_DEFINITIONS
      _CRT_SECURE_NO_WARNINGS
      _SCL_SECURE_NO_WARNINGS
      COMPILED_FROM_DSP
      XML_BUILD_EXPAT
      XML_STATIC
    )
else (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_property(
      DIRECTORY
      PROPERTY
      COMPILE_DEFINITIONS
      _FILE_OFFSET_BITS=64
      _LARGEFILE_SOURCE
      GCC
      HAVE_EXPAT_CONFIG_H
    )
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

include_directories(.)
include_directories(../)
include_directories(../MSToolkit/include)

add_library(
  Comet
  STATIC
  Comet.cpp
  Threading.cpp
  CometSearch.cpp
  CometPreprocess.cpp
  CometMassSpecUtils.cpp
  CometPostAnalysis.cpp
  CometWriteOut.cpp
  CometWriteSqt.cpp
  CometWritePepXML.cpp
)
