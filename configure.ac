AC_INIT(cafe, 4.1)

m4_include([m4/ax_lib_readline.m4])

AC_CONFIG_HEADER([config.h])

AC_LANG(C++)
AC_PROG_CXX([pgc++ icc g++])
AC_PROG_CC([pgcc icc gcc])
AC_OPENMP

ax_blas_ok=no

# Get an appropriate command line for PGI compiler
AC_SUBST(ALL_WARNINGS,-Wall)
wallswitch=`$CC -Wall 2>&1`
if echo "$wallswitch" | grep -i "Error-Unknown" >/dev/null ; then
        AC_SUBST(ALL_WARNINGS,-Minform=inform)
else
        AC_SUBST(C_STANDARD,-std=c11)
fi

AC_SEARCH_LIBS(sgemm, mkl_intel_lp64,
        [ax_blas_ok=yes;BLAS_LIBS="-lmkl_sequential -lmkl_core"],,
        [-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread])

AC_SEARCH_LIBS(cblas_dgemm, satlas,
        [ax_blas_ok=yes;BLAS_LIBS="-lsatlas";AC_DEFINE([HAVE_ATLAS],[1], [Atlas for matrix multiplication])],,
        [])

if test $ax_blas_ok = yes; then
  save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
  AC_DEFINE([HAVE_BLAS],[1], [If a matrix multiplication library exists])
fi

AC_CHECK_HEADERS_ONCE([dirent.h ndir.h])
AX_LIB_READLINE

dnl Process Makefile.in to create Makefile
AC_OUTPUT(Makefile)
