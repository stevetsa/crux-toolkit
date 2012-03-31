# Automake template for hk2ms2 library


AM_CXXFLAGS = -Wall \
              -W \
              $(USEDOUBLES) \
              $(USEDARWIN) \
              -D_LARGEFILE_SOURCE=1 \
              -D_FILE_OFFSET_BITS=64 \
              -D_GNU_SOURCE


libhk2ms2_la_SOURCES = \
	HK2MS2.cpp \
        HKParser.cpp  \
        HKParser.h \
	CruxHK2MS2Application.h \
        CruxHK2MS2Application.cpp

libhk2ms2_la_CXXFLAGS = $(AM_CXXFLAGS) \
  -I$(top_srcdir)/src/c/ \
  -I$(top_srcdir)/src/external/MSToolkit/include \
  -I$(top_srcdir)/src/c/hardklor \
  -DGCC \
  -D_LARGEFILE_SOURCE \
  -DHAVE_EXPAT_CONFIG_H \
  -D_FILE_OFFSET_BITS=64 \
  -DCRUX


noinst_LTLIBRARIES = libhk2ms2.la

