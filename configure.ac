AC_INIT([RWiener], [0.1-0])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS` \ 
        -Wall -pedantic
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`

AC_PROG_CC

AC_CHECK_HEADERS([R], [], [AC_MSG_FAILURE("R header file not found")])
AC_CHECK_LIB([R], [], [AC_MSG_FAILURE("R library not found")])

AC_OUTPUT(src/Makevars)
