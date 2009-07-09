AC_INIT([libmoleculizer], [1.1.2], [nathan.addy@gmail.com])
AM_INIT_AUTOMAKE([-Wall -Werror])

# Ensure this file is present for sanity.
AC_CONFIG_SRCDIR([src/libmoleculizer/mzr/moleculizer.hh]) 


# Make sure we've got compilers and other standard features and such.
# It seems a little wordy though... 

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
AC_FUNC_STRTOD # Checks for library functions.
AC_CHECK_FUNCS([floor modf pow sqrt strtol])


# Here we make sure the mandatory libxml++ is installed.  Because
# libmoleculizer only uses its basic features, we can link in either
# version.

AC_MSG_CHECKING([Searching for libxml++.])
PKG_CHECK_EXISTS( libxml++-2.6, [LINK_WITH_LIBXMLPP_26="yes"], [PKG_CHECK_EXISTS( libxml++-1.0, [LINK_WITH_LIBXMLPP_10="yes"], [NO_LIBXMLPP="yes"])])

if test "$LINK_WITH_LIBXMLPP_26" = "yes"; then
   AC_MSG_RESULT([Found libxml++-2.6.])
   PKG_CHECK_MODULES([LIBXMLPP], libxml++-2.6)
   LIBMZR_REQUIREMENTS="libxml++-2.6"

elif test "$LINK_WITH_LIBXMLPP_10" = "yes"; then
   AC_MSG_RESULT([Found libxml++-1.0.])
   PKG_CHECK_MODULES([LIBXMLPP], libxml++-1.0)
   LIBMZR_REQUIREMENTS="libxml++-1.0"

elif test "$NO_LIBXMLPP" = "yes"; then
    AC_MSG_RESULT([No libxml++ found.])
    AC_MSG_ERROR([Error. You must have some version of libxml++ installed on your system to build libmoleculizer.])
fi

AC_SUBST(LIBXMLPP_CFLAGS)
AC_SUBST(LIBXMLPP_LIBS)
AC_SUBST(LIBMZR_REQUIREMENTS) # For libmoleculizer.pc


## Make sure Python is installed, so that we can embed it into the
## library.
AZ_PYTHON_DEFAULT( )
AZ_PYTHON_PATH( )
AZ_PYTHON_VERSION_ENSURE( [2.2] )		
AZ_PYTHON_CSPEC
AZ_PYTHON_LSPEC

dnl ## Check whether we should build the demo programs.
dnl AC_MSG_CHECKING([whether to build Demos])

dnl AC_ARG_ENABLE([demos],
dnl   AC_HELP_STRING([--build-library-only],
dnl                  [only build the libmoleculizer library, with none of the demo programs]), [disable_demos=yes])

dnl AM_CONDITIONAL([ENABLE_DEMOS], [test ! "$disable_demos" = "yes"])
dnl if test "$disable_demos" = "yes"; then
dnl   AC_MSG_RESULT([no])
dnl else
dnl   AC_MSG_RESULT([yes])
dnl fi

AM_CONDITIONAL([ENABLE_DEMOS], [test "yes" = "yes"])

# Produce this header file.
AC_CONFIG_HEADERS([src/libmoleculizer/moleculizer_config.hh])

AC_CONFIG_FILES([Makefile
                src/Makefile
		src/libmoleculizer/Makefile
		src/libmoleculizer/nauty/Makefile
                src/libmoleculizer/cpx/Makefile
		src/libmoleculizer/dimer/Makefile
		src/libmoleculizer/fnd/Makefile
		src/libmoleculizer/ftr/Makefile
		src/libmoleculizer/mol/Makefile
		src/libmoleculizer/mzr/Makefile
		src/libmoleculizer/nmr/Makefile
		src/libmoleculizer/plex/Makefile
		src/libmoleculizer/stoch/Makefile
		src/libmoleculizer/utl/Makefile
                doc/Makefile
                demos/Makefile
                demos/demo-programs/Makefile
                demos/sample-models/Makefile
		python-src/Makefile
		python-src/language_parser/Makefile
                libmoleculizer-1.0.pc
                doc/user-manual.tex
		])


# Actually calculate and produce the files specified above.
AC_OUTPUT

