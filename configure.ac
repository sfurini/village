# Copyright (C) 2006 Simone Furini <furini.simone@gmail.com>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

######################################################################
#   Initialize                                                       #
######################################################################

AC_INIT([village],
        [4.0],
        [Simone Furini furini.simone@gmail.com],
        [village])
AC_CONFIG_AUX_DIR(config)
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE(village,4.0)

######################################################################
#   Checks for programs                                              #
######################################################################
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB
if test "x$GCC" = "xyes" ; then
	CFLAGS="-g -Wall -O3 -funroll-loops"
	CXXFLAGS="-g -Wall -O3 -funroll-loops"
fi

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])

# Checks for header files.
AC_CHECK_HEADERS([fenv.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL

######################################################################
#   Checks for typedefs, structures, and compiler characteristics    #
######################################################################
AC_C_CONST

######################################################################
#   Check for library functions                                      #
######################################################################
AC_CHECK_FUNCS([floor pow sqrt strchr strrchr feenableexpect])

######################################################################
#   debug support                                                    #
######################################################################
AC_ARG_ENABLE(debug,
	      [  --enable-debug          turn on debugging [default=no]],
	      AC_DEFINE(DEBUG, 1, [Define to 1 to enable debugging code.]),)

######################################################################
#   configure output                                                 #
######################################################################
AC_CONFIG_FILES([
   Makefile
   README
   m4/Makefile
   doc/Makefile
   dat/Makefile
   include/Makefile
   src/Makefile
])
AC_OUTPUT

