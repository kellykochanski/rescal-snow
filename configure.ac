#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.64])
AC_INIT([ReSCAL], [1.6], [rozier@ipgp.fr])
AC_CONFIG_SRCDIR([src])
#AC_CONFIG_SRCDIR([genesis])
#AC_CONFIG_SRCDIR([bin2png])
#AC_CONFIG_SRCDIR([regenesis])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
PKG_CHECK_MODULES([GTK], [gtk+-2.0])
PKG_CHECK_MODULES([GLIB], [glib-2.0])
PKG_CHECK_MODULES([GMODULE], [gmodule-2.0])
PKG_CHECK_MODULES([GTHREAD], [gthread-2.0])
AC_CHECK_LIB([pthread], [pthread_create], [], AC_MSG_ERROR([Required library pthread not found. Install pthread development package and try again.]))
AC_CHECK_LIB([png], [png_write_row], [], AC_MSG_ERROR([Required library libpng not found. Install libpng development package and try again.]))
#AC_CHECK_LIB([gd], [gdImageCreate], [], AC_MSG_ERROR([Required library gd not found. Install gd development package and try again.]))

# Checks for header files.
AC_CHECK_HEADERS([locale.h memory.h stdlib.h string.h sys/time.h unistd.h])

# Checks for functions.
AC_CHECK_FUNCS([strstr])
AC_CHECK_FUNCS([strdup])

# Checks for OpenMP
AC_OPENMP

# Checks for AVX
AC_ARG_ENABLE(avx, [AC_HELP_STRING([--enable-avx],[enable AVX optimizations])], have_avx=$enableval, have_avx=no)
if test "$have_avx" = "yes"; then
  AC_DEFINE(USE_AVX,1,[Define to enable AVX optimizations.])
  AC_SUBST(AVX_CFLAGS,"-mavx")
fi
AM_CONDITIONAL(HAVE_AVX, test "$have_avx" = "yes")

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor gettimeofday memset pow setlocale sqrt strchr])

#AC_CONFIG_FILES([Makefile rescal/Makefile regenesis/Makefile])
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
AC_C_INLINE
AC_PROG_LN_S
