AC_PREREQ(2.59)
AC_INIT([moleculizer], [2.0], [addy@molsci.org])
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

# We should search for and export XML_INCLUDE_DIR and XML_LIB_DIR here.
# As well as anything else.

# We should search for and export * and * here.

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

# Output this header file.
AC_CONFIG_HEADERS([src/moleculizer_config.hh])

AC_OUTPUT

