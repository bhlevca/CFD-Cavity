## process this file with automake to produce Makefile.in

SHELL = /bin/sh
NAVSRCDIR=../../src
DEFS += -DNOPARALLEL

include $(top_srcdir)/src/Makefile.common

bin_PROGRAMS = vrml2nav
vrml2nav_SOURCES = \
	$(NAVHEADERS) \
	$(TOKENIZERSRC) \
	vrml2nav.cxx \
	vrml2nav.h
vrml2nav_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/src
