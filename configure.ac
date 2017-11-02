AC_INIT(cafe, 4.1)
dnl Switch to a C++ compiler, and check if it works.

AC_LANG(C++)
AC_PROG_CXX([icc g++])
AC_PROG_CC([icc gcc])
AC_OPENMP

ax_blas_ok=no

AC_SEARCH_LIBS(sgemm, mkl_intel_lp64,
	[ax_blas_ok=yes;BLAS_LIBS="-lmkl_sequential -lmkl_core"],,
	[-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread])

if test $ax_blas_ok = yes; then
  save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
  AC_DEFINE(HAVE_BLAS,1)
fi

AC_CHECK_HEADERS_ONCE([dirent.h ndir.h])

dnl Process Makefile.in to create Makefile
AC_OUTPUT(Makefile)

