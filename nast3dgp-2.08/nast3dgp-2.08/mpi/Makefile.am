## process this file with automake to produce Makefile.in

@SET_MAKE@
SHELL = /bin/sh
NAVSRCDIR = ../src

# external variables
CXX	= @MPICXX@

if MPIVERSION

include $(top_srcdir)/src/Makefile.common

bin_PROGRAMS = navcalcmpi
navcalcmpi_SOURCES = $(NAVHEADERS) $(COMMONSRC) $(CALCSRC)
navcalcmpi_LDADD = @MPILIBS@

endif
