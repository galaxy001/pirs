AC_PREREQ([2.59])
AC_INIT([pirs], [2.0.0])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([pirs/main.cpp])
AC_CONFIG_HEADER([config.h])

AC_LANG([C++])
AC_PROG_CXX
AC_PROG_CC

AC_CHECK_LIB([z], [gzopen])

AC_MSG_CHECKING([Whether to include support for 'pirs diploid' command])
AC_ARG_ENABLE([pirs-diploid],
	AS_HELP_STRING([--disable-pirs-diploid], [Do not support the 'pirs diploid' command]),
	[enable_pirs_diploid=$enableval],
	[enable_pirs_diploid=yes])
if test "$enable_pirs_diploid" = "yes"; then
	AC_DEFINE([ENABLE_PIRS_DIPLOID], [1], [Define to 1 to include support for
		   the 'pirs diploid' command])
fi
AM_CONDITIONAL([ENABLE_PIRS_DIPLOID], [test "$enable_pirs_diploid" = "yes"])
AC_MSG_RESULT([$enable_pirs_diploid])

AC_MSG_CHECKING([Whether to include support for 'pirs simulate' command])
AC_ARG_ENABLE([pirs-simulate],
	AS_HELP_STRING([--disable-pirs-simulate], [Do not support the 'pirs simulate' command]),
	[enable_pirs_simulate=$enableval],
	[enable_pirs_simulate=yes])
if test "$enable_pirs_simulate" = "yes" ; then
	AC_DEFINE([ENABLE_PIRS_SIMULATE], [1], [Define to 1 to include support for
		   the 'pirs simulate' command])
fi
AM_CONDITIONAL([ENABLE_PIRS_SIMULATE], [test "$enable_pirs_simulate" = "yes"])
AC_MSG_RESULT([$enable_pirs_simulate])

AC_MSG_CHECKING([Whether to include support for multiple threads])
AC_ARG_ENABLE([threads],
	AS_HELP_STRING([--disable-threads],
		       [Build a single-threaded program]),
	[ENABLE_THREADS=$enableval],
	[ENABLE_THREADS=yes])
if test $ENABLE_THREADS = yes; then
	AC_DEFINE([ENABLE_THREADS], [1],
		  [Define to 1 to include support for multiple threads])
	AC_MSG_RESULT([$ENABLE_THREADS])
	AC_OPENMP
else
	AC_MSG_RESULT([$ENABLE_THREADS])
fi
AM_CONDITIONAL([ENABLE_THREADS], [test $ENABLE_THREADS = yes])

AC_SUBST([OPENMP_LDFLAGS], [$OPENMP_LDFLAGS])
AC_SUBST([OPENMP_CXXFLAGS], [$OPENMP_CXXFLAGS])

AC_MSG_CHECKING([Whether to include support for SSE2 instructions])
AC_ARG_ENABLE([sse2],
	AS_HELP_STRING([--disable-ssse2],
		       [Use non-SSE2 version of Mersenne Twister]),
	[HAVE_SSE2=$enableval],
	[HAVE_SSE2=yes])
if test $HAVE_SSE2 = yes; then
	AC_DEFINE([HAVE_SSE2], [1],
		  [Define to 1 to include support for SSE2 instructions])
	AC_SUBST([RNG_CPPFLAGS], ["-DSFMT_MEXP=19937"])
	AC_SUBST([RNG_CXXFLAGS], ["-msse2"])
fi
AM_CONDITIONAL([HAVE_SSE2], [test $HAVE_SSE2 = yes])
AC_MSG_RESULT([$HAVE_SSE2])

AC_OUTPUT([
Makefile
pirs/Makefile
])

