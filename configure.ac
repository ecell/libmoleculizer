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
dnl
dnl PKG_CHECK_MODULES(LIBXML, libxml++-1.0)
PKG_CHECK_MODULES(LIBXML, libxml++-2.6 >= 2.10)
AC_SUBST(LIBXML_CFLAGS)
AC_SUBST(LIBXML_LIBS)

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


AC_ARG_ENABLE(static-build,
 [  --enable-static-build  Build the entire project as one monolithic staticly linked program],
 [ static_build=yes
 ])

AM_CONDITIONAL([STATIC_BUILD], [test "$static_build" = yes])

AC_ARG_ENABLE(index-alphabetically,
        [ --enable-alphabetic-indexing Causes binding sites in the naming unit to be indexed alphabetically.]
        [index_alphabetically=yes])



# Output this header file.
AC_CONFIG_HEADERS([src/moleculizer_config.hh])

AC_OUTPUT

