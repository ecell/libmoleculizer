AC_INIT([libmoleculizer], [1.1.4], [nathan.addy@gmail.com])
AM_INIT_AUTOMAKE([-Wall -Werror])

INSTALLATION_VERSION=1.0
AC_SUBST(INSTALLATION_VERSION)

# Ensure this file is present for sanity.
AC_CONFIG_SRCDIR([src/libmoleculizer/mzr/moleculizer.hpp]) 

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

# Make sure the libxml++ is installed.  We look for, and can link in either
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
AZ_PYTHON_PATH( )
AZ_PYTHON_VERSION_ENSURE( [2.3] )		
AZ_PYTHON_CSPEC
AZ_PYTHON_LSPEC

if test ! "$prefix" = "NONE"; then
   echo "Installing python to directory $prefix"
   PYTHON_SETUP_PY_PREFIX_FLAGS="-f --prefix=$prefix"
else
   echo "Installing python to al location"
   PYTHON_SETUP_PY_PREFIX_FLAGS= 
fi

AC_SUBST(PYTHON_SETUP_PY_PREFIX_FLAGS)

## Check whether we should build the demo programs.
AC_MSG_CHECKING([whether to build Demos])

AC_ARG_ENABLE([demos],
   AC_HELP_STRING([--disable-demos],
                  [only build the libmoleculizer library, with none of the demo programs]), [demos=no], [demos=yes])

AM_CONDITIONAL([ENABLE_DEMOS], [test "$demos" = "yes"])
if test "$demos" = "yes"; then
   AC_MSG_RESULT([yes])
 else
   AC_MSG_RESULT([no])
 fi

# Produce this header file.
AC_CONFIG_HEADERS([src/libmoleculizer/moleculizer_config.hpp])

AC_CONFIG_FILES([Makefile
                documentation/Makefile
                libmoleculizer-1.0.pc
                documentation/user-manual.tex
                src/bindings/python/Makefile
                src/demo-programs/Makefile
                src/Makefile
                src/libmoleculizer/cpx/Makefile
                src/libmoleculizer/plex/Makefile
                src/libmoleculizer/stoch/Makefile
                src/libmoleculizer/nauty/Makefile
                src/libmoleculizer/Makefile
                src/libmoleculizer/ftr/Makefile
                src/libmoleculizer/dimer/Makefile
                src/libmoleculizer/mol/Makefile
                src/libmoleculizer/utl/Makefile
                src/libmoleculizer/template/tests/Makefile
                src/libmoleculizer/fnd/Makefile
                src/libmoleculizer/nmr/Makefile
                src/libmoleculizer/nmr/tests/Makefile
                src/libmoleculizer/mzr/Makefile
                src/libmoleculizer/mzr/tests/Makefile
                src/src-python/Makefile
                src/src-python/libmoleculizer/mzr_parser/Makefile
                examples/Makefile
                ])



# Actually calculate and produce the files specified above.
AC_OUTPUT

