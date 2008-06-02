AC_PREREQ(2.59)
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

dnl PKG_CHECK_MODULES(LIBXML, libxml++-1.0)
PKG_CHECK_MODULES(LIBXML, libxml++-2.6)
AC_SUBST(LIBXML_CFLAGS)
AC_SUBST(LIBXML_LIBS)

PKG_CHECK_MODULES(GSL, gsl)
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

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
		src/stoch/Makefile
		src/utl/Makefile
		])


AC_ARG_ENABLE(static-build,
 [  --enable-static-build  Build the entire project as one monolithic staticly linked program],
 [ static_build=yes
 ])

AM_CONDITIONAL([STATIC_BUILD], [test "$static_build" = yes])

# Output this header file.
AC_CONFIG_HEADERS([src/moleculizer_config.hh])

AC_OUTPUT

