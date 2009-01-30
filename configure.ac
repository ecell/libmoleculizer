AC_INIT([libmoleculizer], [1.0], [addy@molsci.org])
AM_INIT_AUTOMAKE([-Wall -Werror])

# Ensure this file is present for sanity.
AC_CONFIG_SRCDIR([src/mzr/moleculizer.hh]) 

# This will be added at alpha-version.
# AC_CONFIG_TESTDIR([tests])

# This should be be greatly improved.
ACLOCAL="aclocal"
AC_SUBST( ACLOCAL )

# Do I need to add somethng about dynamic linking?
AC_PROG_CC			
AC_PROG_CXX
AC_PROG_CPP
AC_PROG_LIBTOOL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_INSTALL

AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h])

AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor modf pow sqrt strtol])

dnl
dnl
dnl
AM_PATH_PYTHON(2.3)
AC_CHECK_PROGS(PYCHECKER, [pychecker], [no])
if test "$PYCHECKER" == no ; then
   AC_MSG_WARN([pychecker not found.]) 
fi

dnl
dnl
AC_CHECK_LIB(python${PYTHON_VERSION},main,,AC_MSG_ERROR([could not find libpython.]))
dnl AC_CHECK_LIB(boost_python,main,,AC_MSG_ERROR([could not find libboost_python.]),-lpython${PYTHON_VERSION})
dnl

PKG_CHECK_EXISTS( libxml++-2.6, [PKG_CHECK_MODULES(LIBXML, libxml++-2.6 >= 2.10)], [PKG_CHECK_MODULES(LIBXML, libxml++-1.0)])
PKG_CHECK_EXISTS( libxml++-2.6, [LIBXMLPP_VERSION="2.6"], [LIBXMLPP_VERSION="1.0"])

AC_SUBST(LIBXML_CFLAGS)
AC_SUBST(LIBXML_LIBS)
AC_SUBST(LIBXMLPP_VERSION)

# Checks for libraries.
AC_CONFIG_FILES([Makefile
                src/Makefile
                src/cpx/Makefile
		src/dimer/Makefile
		src/fnd/Makefile
		src/ftr/Makefile
		src/mol/Makefile
		src/mzr/Makefile
		src/nmr/Makefile
		src/plex/Makefile
                src/python/Makefile
		src/stoch/Makefile
		src/utl/Makefile
                src/nmr/tests/Makefile
                src/mzr/tests/Makefile
                doc/Makefile
                doc/manual/Makefile
                doc/demos/Makefile
                doc/demos/simulators/Makefile
                libmoleculizer-1.0.pc
		])


AC_ARG_ENABLE([demos],
  AC_HELP_STRING([--disable-demos],
                 [build libmoleculizer without the associated demos.]), [],
  [enable_demos=yes]
)

AC_MSG_CHECKING([whether to build Demos])
if test "$enable_demos" = "yes"; then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

AM_CONDITIONAL([ENABLE_DEMOS], [test "$enable_demos" = "yes"])

# Output this header file.
AC_CONFIG_HEADERS([src/moleculizer_config.hh])

AC_OUTPUT

