#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT([NaSt3DGP], [2.08], [nast3dgp@ins.uni-bonn.de])
AC_CONFIG_SRCDIR([src/navier.h])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_HEADER([src/config.h])
AM_INIT_AUTOMAKE([1.7 -Wall foreign])
AM_MAINTAINER_MODE

AC_LANG([C++])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
dnl FIXME: Fortran for ./util/GridGen/GridGen.f
AC_PROG_MAKE_SET
AC_PROG_LN_S

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
dnl AC_FUNC_MALLOC
dnl AC_FUNC_REALLOC
dnl AC_FUNC_STRTOD
dnl AC_CHECK_FUNCS([sqrt strchr strtol])

# Checks for system services
AC_PATH_XTRA

# Save the cache since the following tests might fail more often
AC_CACHE_SAVE

# TIMED
AC_ARG_ENABLE([timer],
	      AC_HELP_STRING([--enable-timer],
			     [output timing information (no)]),
	      [enable_timer=$enableval], [enable_timer=no])
if test x"$enable_timer" = xyes
then
	AC_DEFINE([TIMED], 1, [Define if timing information is wanted.])
fi

#DEBUG
#notice: when enabling debug, compile flags are set equal to -g
# that is, user given flags are overridden. This is not nice.
AC_ARG_ENABLE([debug],
		AC_HELP_STRING([--enable-debug],
				[compile with debug information (no)]),
				[enable_debug=$enableval], [enable_debug=no])
if test x"$enable_debug" = xno
then
	AC_MSG_NOTICE([debugging disabled])
	AC_DEFINE([NDEBUG], 1, [Define to disable debugging.])
else
	AC_MSG_NOTICE([enabling debug info])
	CFLAGS="-g"
	CXXFLAGS="-g"
fi

# MPI
AC_ARG_ENABLE([mpi],
	      AC_HELP_STRING([--enable-mpi],
			     [build MPI parallel version
			      (default is yes if MPI code can be compiled)]),
	      [enable_mpi=$enableval], [enable_mpi=maybe])
disabled_msg="MPI version disabled."
if test x"$enable_mpi" = xyes -o x"$enable_mpi" = xmaybe
then
	AC_MSG_NOTICE([checking for MPI features])
	ACX_MPI([enable_mpi=yes],
		[err_msg="cannot find out how to compile MPI code."
		 if test x"$enable_mpi" = xyes
		 then
			 AC_MSG_ERROR([$err_msg])
		 else
			 AC_MSG_WARN([$err_msg $disabled_msg])
			 enable_mpi=no
		 fi])
else
	AC_MSG_NOTICE([$disabled_msg])
fi
AC_ARG_VAR(MPILIBS,[necessary libraries to link MPI C code])
AM_CONDITIONAL(MPIVERSION, [test x"$enable_mpi" = xyes])

# VTK
# FIXME: need to check for presence of VTK headers and libraries.
#AC_ARG_ENABLE([vtk],
#	    AC_HELP_STRING([--enable-vtk@<:@=VTK_ROOT@:>@],
#			   [enable VTK visualization tool.
#			    Argument VTK_ROOT can be supplied as the installation
#			    base directory of the VTK headers and libraries
#			    if they are not found in standard locations.
#			    Default is no.]),
#	    [enable_vtk=$enableval], [enable_vtk=no])
#AC_MSG_NOTICE([enable_vtk is $enable_vtk])
#if test x"$enable_vtk" != xno -a x"$enable_vtk" != xyes
#then
#	VTK_ROOT="$enable_vtk"
#else
	# FIXME
#	VTK_ROOT=/usr
#fi
#AC_SUBST(VTK_ROOT)
#AM_CONDITIONAL(VTK, [test x"$enable_vtk" != xno])


dnl Unfortunately, autoconf lacks a --docdir option. This is a workaround:
AC_ARG_WITH([docdir],
	    AC_HELP_STRING([--with-docdir=DIR],
			   [Where to put the documentation @<:@PREFIX/share/doc/$PACKAGE@:>@]),
	    [if test "x$withval" = xyes || test "x$withval" = xno; then
	       docdir='${prefix}/share/doc/'
	     else
               docdir="$withval"
	     fi],
	     [docdir='${prefix}/share/doc/'"$PACKAGE"])
AC_SUBST(docdir)

AC_CONFIG_FILES([Makefile
                 src/Makefile
		 doc/Makefile
		 tools/Makefile
		 tools/vrml2nav/Makefile
		 mpi/Makefile])
AC_OUTPUT

