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
AC_PROG_LIBTOOL
AC_PROG_MAKE_SET
AC_PROG_INSTALL

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
AC_CHECK_LIB(boost_python,main,,AC_MSG_ERROR([could not find libboost_python.]),-lpython${PYTHON_VERSION})
dnl
dnl
dnl PKG_CHECK_MODULES(LIBXML, libxml++-1.0)
PKG_CHECK_MODULES(LIBXML, libxml++-2.6)
AC_SUBST(LIBXML_CFLAGS)
AC_SUBST(LIBXML_LIBS)

# Checks for libraries.
AC_CONFIG_FILES([Makefile
                src/Makefile
		src/nauty/Makefile
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
                doc/demos/simulators/particle/src/Makefile
                doc/demos/simulators/stochastic/src/Makefile
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

