AC_INIT([libmoleculizer], [1.1.1], [nathan.addy@gmail.com])
AM_INIT_AUTOMAKE([-Wall -Werror])

# Ensure this file is present for sanity.
AC_CONFIG_SRCDIR([src/libmoleculizer/mzr/moleculizer.hh]) 

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

AC_MSG_CHECKING([Searching for libxml++.])
PKG_CHECK_EXISTS( libxml++-2.6, [LINK_WITH_LIBXMLPP_26="yes"], [PKG_CHECK_EXISTS( libxml++-1.0, [LINK_WITH_LIBXMLPP_10="yes"], [NO_LIBXMLPP="yes"])])
if test "$LINK_WITH_LIBXMLPP_26" = "yes"; then
   AC_MSG_RESULT([Found 2.6.])
   PKG_CHECK_MODULES([LIBXMLPP], libxml++-2.6)
   LIBMZR_REQUIREMENTS="libxml++-2.6"

elif test "$LINK_WITH_LIBXMLPP_10" = "yes"; then
   AC_MSG_RESULT([Found 1.0.]
   PKG_CHECK_MODULES([LIBXMLPP], libxml++-1.0)
   LIBMZR_REQUIREMENTS="libxml++-1.0"
else
    AC_MSG_RESULT([None found.])
    AC_MSG_ERROR([Error. You must have libxml++ installed on your system to build libmoleculizer.])
fi

AC_SUBST(LIBXMLPP_CFLAGS)
AC_SUBST(LIBXMLPP_LIBS)
AC_SUBST(LIBMZR_REQUIREMENTS)

## Search for Python and add it to everything.
AZ_PYTHON_DEFAULT
AZ_PYTHON_CSPEC
AZ_PYTHON_LSPEC

AC_SUBST(PYTHON_CPPFLAGS)
AC_SUBST(PYTHON_LDFLAGS)
AC_SUBST(PYTHON_EXTRA_LIBS)
AC_SUBST(PYTHON_EXTRA_LDFLAGS)

## Determine whether or not to build the demo programs
AC_ARG_ENABLE([demos],
  AC_HELP_STRING([--library-only],
                 [only build libmoleculizer, and none of the extra demo programs]), [disable_demos=yes])

AC_MSG_CHECKING([whether to build Demos])
AM_CONDITIONAL([ENABLE_DEMOS], [test ! "$disable_demos" = "yes"])
if test "$disable_demos" = "yes"; then
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([yes])
fi

# Output this header file.
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

AC_OUTPUT

