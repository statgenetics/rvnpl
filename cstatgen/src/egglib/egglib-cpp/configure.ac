#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([EggLib C++ library], [2.1.7], [demita@gmail.com], [egglib-cpp])
: ${CFLAGS=""}  # no -g -O2 option
: ${CXXFLAGS=""}  # no -g -O2 option
AC_PREREQ([2.63])
AC_CONFIG_SRCDIR([SitePolymorphism.hpp])
AC_CONFIG_HEADERS([config.h])
AC_PROG_RANLIB

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB(bpp-core,    main, , )
AC_CHECK_LIB(bpp-seq,     main, , )
AC_CHECK_LIB(bpp-popgen,  main, [bpp_available=true], [bpp_available=false])
AM_CONDITIONAL([BPP_AVAILABLE], [test x$bpp_available = xtrue])
AC_CHECK_LIB(gslcblas,  main)
AC_CHECK_LIB(gsl,  main, [gsl_available=true], [gsl_available=false])
AM_CONDITIONAL([GSL_AVAILABLE], [test x$gsl_available = xtrue])
#AC_OPENMP

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
#AC_FUNC_MALLOC
#AC_FUNC_REALLOC
AC_CHECK_FUNCS([sqrt])
AC_TYPE_SIZE_T

AM_INIT_AUTOMAKE
AC_CONFIG_FILES(Makefile)
AC_OUTPUT

