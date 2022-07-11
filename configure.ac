#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([H2OI95], [1.1], [cory.marshall@drdc-rddc.gc.ca])
AM_INIT_AUTOMAKE
LT_PREREQ([2.4.7])
LT_INIT
AC_CONFIG_SRCDIR([src/H2OI95.f])
# Checks for programs.
AC_PROG_F77
AC_PROG_INSTALL
AC_PROG_SED
AC_CHECK_PROG([DIFF],[diff],[yes])
PKG_INSTALLDIR

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_TESTDIR([tests])
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/libH2OI95/Makefile
                 src/libH2OI95/libH2OI95.pc
                 tests/Makefile])
AC_OUTPUT
