AC_INIT
AC_PREREQ(2.59)
AC_CONFIG_SRCDIR(src/c/ion_series.c)
AM_CONFIG_HEADER([config.h])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE(crux, 1.0.0)
AC_PREFIX_DEFAULT(${HOME})

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PATH_PROG(LD, ld)
AC_PATH_PROG(SED, sed)
AC_PATH_PROG(PYTHON, python)

AC_DEFINE_UNQUOTED(HOSTNAME, "`hostname`", [Build host name])

# Checks for libraries.
AC_CHECK_LIB([m], 
             [isnan, sqrt, log, exp, fabs, pow, floor, ceil], 
             [], 
             [AC_MSG_WARN([Math library not found])])

# Checks for header files.
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS([arpa/inet.h ctype.h fcntl.h limits.h netdb.h netinet/in.h stdlib.h string.h strings.h sys/socket.h sys/time.h unistd.h ])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_HEADER_STDBOOL
AC_TYPE_PID_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_TYPE_SIGNAL
AC_FUNC_STAT
AC_CHECK_FUNCS([floor getcwd gethostbyname isascii memset pow rint socket sqrt strchr strcspn strdup strspn strstr])
AC_CHECK_FUNCS([strlcpy], [have_strlcpy=yes], [have_strlcpy=no])
AM_CONDITIONAL(WANT_STRLCPY, test ["$have_strlcpy" = no])

# who is installing 
USER=`whoami`

# Check for platform-specific settings
case $host in
   *cygwin*)
      AC_DEFINE(Cygwin, [], [Using cygwin])
      MACHID="cygwin"
      ;;
   *darwin*)
      AC_DEFINE(MACOSX, [], [Using MacOS X])
      MACHID="darwin"
	  have_startup="yes"
      ;;
   *linux*)
      AC_DEFINE(Linux, [], [Using Linux])
      MACHID="linux"
	  have_startup="yes"
      ;;
   *aix*)
      ac_mpicc_flags="-O3 -qstrict -qarch=pwr3 -qtune=pwr3"
      AC_DEFINE(AIX, [], [Using AIX])
      MACHID="aix"
      if /usr/sbin/lsattr -EHl proc0 | grep POWER >/dev/null 2>&1; then
          MACHID=rs6000
      fi
      ;;
   *cray*)
      ac_mpicc_flags="-O3 -h msglevel_3"
      AC_DEFINE(cray3e, [], [Using cray3e])
      MACHID="cray"
      ;;
   *irix6*)
      ac_mpicc_flags="-O3"
      AC_DEFINE(IRIX64, [], [Using IRIX64])
      MACHID="irix"
	  have_startup="yes"
      ;;
   *irix*)
      ac_mpicc_flags="-O3"
      MACHID="irix"
	  have_startup="yes"
      ;;
   *solaris*|*sunos*)
      if test "${CC}" != "gcc"; then
          ac_cc_opt="-O -Xa"
      fi
      AC_DEFINE(SunOS, [], [Using SunOS])
      AC_CHECK_LIB([nsl], [gethostbyname], [], [AC_MSG_ERROR([Nsl library not found])])
      AC_CHECK_LIB([socket], [socket], [], [AC_MSG_ERROR([Socket library not found])])
      LIBS="-lsocket -lnsl -lm"
      MACHID="solaris"
	  have_startup="yes"
      ;;
esac

AM_CONDITIONAL(WANT_STARTUP, test ["$have_startup" = yes])
AM_CONDITIONAL(WANT_DARWIN, test ["$MACHID" = darwin])

####################################
# Check for configuration switches
####################################

# debugging
AC_ARG_ENABLE(
    debug,
    [  --enable-debug           Produce an executable with debugging symbols. Default is disabled.],
    [if test "$enableval" = "yes"; then
        DEBUG_FLAGS="-g";
     fi],
    [DEBUG_FLAGS=""])

# optimization
AC_ARG_ENABLE(
    opt,
    [  --enable-opt             Produce an executable with optimization. Default compiler dependent.],
    [if test "$enableval" = "yes"; then
       OPT_FLAGS="$ac_cc_opt";
     fi],
    [OPT_FLAGS="$ac_cc_opt"])

# set percolator location
AC_ARG_WITH(
    percolator-dir,
    [  --with-percolator-dir=DIR  Set path to percolator source. \
      Default: ${top_srcdir}/percolator],
    [PERCOLATOR_DIR=$withval],
    [PERCOLATOR_DIR="\${top_srcdir}/percolator"])

AC_SUBST(MACHID)
AC_SUBST(DEBUG_FLAGS)
AC_SUBST(OPT_FLAGS)
AC_SUBST(LIBS)
AC_SUBST(VERSION)
AC_SUBST(LD)
AC_SUBST(LAMBIN)
AC_SUBST(HOSTNAME)
AC_SUBST(PERCOLATOR_DIR)

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
  Linker:                     ${LD}
  Special Libs:               ${LIBS}
  Python:                     ${PYTHON}
  percolator directory:       ${PERCOLATOR_DIR}

  Run the following commands to compile, test and install meme:
        make   
        make test
        make install
"

