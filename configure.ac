AC_PREREQ([2.68])
AC_INIT([VC2REFERENCE],[0.1.0.2])

AC_CANONICAL_TARGET
case $target in
x86_64-*)
        TARGET="x86_64"
        ;;
*)
        TARGET="unknown"
        ;;
esac
if test "$TARGET" != "x86_64" ; then
   AC_MSG_ERROR([Only x86_64 architecture is supported])
fi;

AS_NANO(VC2REFERENCE_CVS=no,VC2REFERENCE_CVS=yes)

AM_INIT_AUTOMAKE([1.11 silent-rules subdir-objects])
AM_SILENT_RULES([yes])

# Version number that will be appended to the library filename
VC2REFERENCE_MAJORMINOR=0.1
AC_SUBST(VC2REFERENCE_MAJORMINOR)

AC_CONFIG_HEADERS([config.h])

# CURRENT, REVISION, AGE
# - library source changed but interface the same -> increment REVISION
# - interfaces added/removed/changed -> increment CURRENT, set REVISION = 0
# - interfaces added -> increment AGE also
# - interfaces removed -> set AGE = 0
#    (AGE is the number of previous values of CURRENT that are compatible)
VC2REFERENCE_LIBVERSION="0:0:0"
AC_SUBST(VC2REFERENCE_LIBVERSION)

LT_PREREQ([2.2.6])
LT_INIT(disable-static win32-dll)

AC_PROG_CXX
AM_PROG_CC_C_O

AC_CONFIG_SRCDIR([src/Library/Utils.h])

AC_CONFIG_MACRO_DIR([m4])
ACLOCAL_AMFLAGS="-I m4 $ACLOCAL_AMFLAGS"
AC_SUBST(ACLOCAL_AMFLAGS)

AC_MSG_CHECKING(whether to enable other decoders)
AC_ARG_ENABLE([other-decoders],
[  --enable-other-decoders    Turn on other decoders],
[case "${enableval}" in
  yes) otherdecoders=true ;;
  no)  otherdecoders=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-other-decoders]) ;;
esac],[otherdecoders=false])
AM_CONDITIONAL([ENABLE_OTHER_DECODERS], [test "x$otherdecoders" = "xtrue"])

VC2REFERENCE_CFLAGS="$VC2REFERENCE_CFLAGS -g -Og -DDEBUG"
CXXFLAGS="-Og"
VC2REFERENCE_LDFLAGS="$VC2REFERENCE_LDFLAGS"

LIBBOOST_VER="1.49.0"
AX_BOOST_BASE($LIBBOOST_VER, HAVE_BOOST=yes, HAVE_BOOST=no)
if test "x${HAVE_BOOST}" != xyes ; then
  AC_MSG_ERROR([Boost library is required])
fi
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_THREAD
AX_BOOST_SYSTEM
AX_CXX_HAVE_SSTREAM

# Check for pkg-config
PKG_PROG_PKG_CONFIG([0.26])

IS_LINUX=no
AC_CANONICAL_HOST
case $host_os in
  linux*)
    IS_LINUX=yes
    ;;
esac
AM_CONDITIONAL(IS_LINUX, test x${IS_LINUX} = xyes)
AC_SUBST(IS_LINUX)

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AX_CREATE_STDINT_H([src/vc2reference-stdint.h])

AC_APPEND_SUPPORTED_CFLAGS(VC2REFERENCE_CFLAGS, [-Wall])

VC2REFERENCE_CFLAGS="$VC2REFERENCE_CFLAGS -I\$(top_srcdir)"
AC_SUBST(VC2REFERENCE_CFLAGS)

VC2REFERENCE_LDFLAGS="$VC2REFERENCE_LDFLAGS"
AC_APPEND_SUPPORTED_LDFLAGS(VC2REFERENCE_LDFLAGS, [-Wl,--no-undefined])
if test "x$IS_LINUX" = "xno"
then
  AC_APPEND_SUPPORTED_LDFLAGS(VC2REFERENCE_LDFLAGS, [-Qunused-arguments])
fi
AC_SUBST(VC2REFERENCE_LDFLAGS)

VC2REFERENCE_LIBS="\$(top_builddir)/src/Libaray/libVC2-$VC2REFERENCE_MAJORMINOR.la"
AC_SUBST(VC2REFERENCE_LIBS)

AC_SUBST(VC2REFERENCE_PKG_DEPS)

pkgconfigdir="\$(libdir)/pkgconfig"
AC_SUBST(pkgconfigdir)

vc2referencebindir="\$(libdir)/vc2reference/"
AC_SUBST(vc2referencebindir)

configfiledir="\$(sysconfdir)/vc2reference"
AC_SUBST(configfiledir)

AC_CONFIG_FILES([
Makefile
tools/Makefile
src/Makefile
src/boost/Makefile
src/tclap/Makefile
src/Library/Makefile
src/DecodeStream/Makefile
src/DecodeHQ/Makefile
src/DecodeLD/Makefile
src/EncodeHQ-CBR/Makefile
src/EncodeHQ-ConstQ/Makefile
src/EncodeLD/Makefile
])
AC_OUTPUT
