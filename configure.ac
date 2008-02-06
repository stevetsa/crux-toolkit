# automake/autoconf startup stuff
AC_INIT
AC_CONFIG_SRCDIR(src/c/ion_series.c)
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE(crux, 1.0.0)
AM_CONFIG_HEADER([config.h])

# Check for required programs.
AC_PROG_CC
AC_PROG_LIBTOOL
AC_PATH_PROG(SED, sed)
AC_PATH_PROG(PYTHON, python)

# Checks for required header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([ \
    assert.h \
    ctype.h \
    dirent.h \
    ieeefp.h \
    limits.h \
    stdio.h \
    stdlib.h \
    string.h \
    sys/stat.h \
    time.h \
    unistd.h \
])

# Check for required library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

# Checks for required math functions.
AC_CHECK_LIB(
  [m], 
  [sqrt, log, exp, fabs, pow],
  [], 
  [AC_MSG_WARN([Math library not found])]
)

# Set platform-specific settings
case $host in
   *cygwin*)
      AC_DEFINE(Cygwin, [], [Using cygwin])
      MACHID="cygwin"
      ;;
   *darwin*)
      AC_DEFINE(MACOSX, [], [Using MacOS X])
      MACHID="darwin"
      ;;
   *linux*)
      AC_DEFINE(Linux, [], [Using Linux])
      MACHID="linux"
      ;;
esac

####################################
# Check for configuration switches
####################################

# set percolator location
AC_ARG_WITH(
    percolator-dir,
    [  --with-percolator-dir=DIR  Set path to percolator source. \
      Default: ${top_srcdir}/percolator],
    [PERCOLATOR_DIR=$withval],
    [PERCOLATOR_DIR="\${top_srcdir}/../percolator"])

# These variables will be available in the Makefiles.
AC_SUBST(MACHID)
AC_SUBST(VERSION)
AC_SUBST(PERCOLATOR_DIR)

# List all the Makefiles for the project.
AC_CONFIG_FILES([\
    Makefile \
    doc/Makefile \
    doc/user/Makefile \
    doc/developer/Makefile \
    src/Makefile \
    src/c/Makefile \
])

AC_OUTPUT

echo \
"
        ========================
        Configuration parameters
        ========================

  Install path:               ${prefix}
  Version:                    ${VERSION}
  C compiler:                 ${CC}
  C compiler flags:           ${CFLAGS}
  Python:                     ${PYTHON}
  percolator directory:       ${PERCOLATOR_DIR}

  Run the following commands to compile, test and install meme:
        make   
        make test
        make install
"

