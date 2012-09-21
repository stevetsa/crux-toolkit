include(ExternalProject)

ExternalProject_Add(
  ProteoWizard
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  URL file://${CMAKE_CURRENT_SOURCE_DIR}/proteowizard/pwiz-src-without-tv-3_0_3533.tar.bz2
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ./quickbuild.sh --without-mz5 --prefix=${CMAKE_CURRENT_SOURCE_DIR} libraries pwiz/data/msdata/ pwiz/data/identdata/
  INSTALL_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/proteowizard/create_links.sh --prefix=${CMAKE_CURRENT_SOURCE_DIR}
  BUILD_IN_SOURCE 1
)

