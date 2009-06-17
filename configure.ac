AC_INIT([libmoleculizer], [1.1.0], [addy@molsci.org])
AM_INIT_AUTOMAKE([-Wall -Werror])

# Ensure this file is present for sanity.
AC_CONFIG_SRCDIR([src/libmoleculizer/mzr/moleculizer.hh]) 

# This will be added at alpha-version.
# AC_CONFIG_TESTDIR([tests])

# This should be be greatly improved.
#ACLOCAL="aclocal"
#AC_SUBST( ACLOCAL )

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

PKG_CHECK_EXISTS( libxml++-2.6, [LINK_WITH_LIBXMLPP_26="yes"], [PKG_CHECK_EXISTS( libxml++-1.0, [LINK_WITH_LIBXMLPP_10="yes"], [NO_LIBXMLPP="yes"])])



if test "$LINK_WITH_LIBXMLPP_26" = "yes"; then
   echo "Found libxml++-2.6 installed, using..."
   PKG_CHECK_MODULES([LIBXMLPP], libxml++-2.6)
   LIBMZR_REQUIREMENTS="libxml++-2.6"

elif test "$LINK_WITH_LIBXMLPP_10" = "yes"; then

   echo "Found libxml++-1.0 installed, using..."
   PKG_CHECK_MODULES([LIBXMLPP], libxml++-1.0)
   LIBMZR_REQUIREMENTS="libxml++-1.0"

elif test "$NO_LIBXMLPP" = "yes"; then
     echo "No libxml++ was found, building..."
     build_libxmlpp=yes
     LIBXMLPP_CFLAGS=-I../../../libxml++-1.0.5/
     LIBXMLPP_LIBS=
     LIBMZR_REQUIREMENTS=
     EXTRA_DIRS="libxml++-1.0.5"
     UTL_FLAGS="-I../../../libxml++-1.0.5"
     UTL_LIBS="../../../libxml++-1.0.5/libxml++/libxmlpp.la"
     # Should move this to NO_LIBXMLPP if possible.
     AC_CONFIG_SUBDIRS([libxml++-1.0.5])	

fi

AC_SUBST(LIBXMLPP_CFLAGS)
AC_SUBST(LIBXMLPP_LIBS)
AC_SUBST(LIBMZR_REQUIREMENTS)
AC_SUBST(UTL_FLAGS)
AC_SUBST(UTL_LIBS)
AC_SUBST(EXTRA_DIRS)

# Checks for libraries.
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


AC_ARG_ENABLE([demos],
  AC_HELP_STRING([--enable-demos],
                 [build demo programs along with libmoleculizer]), [enable_demos=yes])

AC_MSG_CHECKING([whether to build Demos])
if test "$enable_demos" = "yes"; then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

AC_PYTHON_DEVEL

AC_SUBST(PYTHON_CPPFLAGS)
AC_SUBST(PYTHON_LDFLAGS)
AC_SUBST(PYTHON_EXTRA_LIBS)
AC_SUBST(PYTHON_EXTRA_LDFLAGS)

AM_CONDITIONAL([ENABLE_DEMOS], [test "$enable_demos" = "yes"])
AM_CONDITIONAL([BUILD_LIBXMLPP], [test "build_libxmlpp" = "yes"])

# Output this header file.
AC_CONFIG_HEADERS([src/libmoleculizer/moleculizer_config.hh])

AC_OUTPUT

