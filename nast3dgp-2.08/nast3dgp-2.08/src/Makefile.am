## process this file with automake to produce Makefile.in

# This Makefile builds the serial components.

SHELL = /bin/sh
NAVSRCDIR = .
DEFS += -DNOPARALLEL

EXTRA_DIST = Makefile.VC5

include $(top_srcdir)/src/Makefile.common

SETUPSRC = \
	$(NAVSRCDIR)/object.cxx \
	$(NAVSRCDIR)/object.h \
	$(NAVSRCDIR)/parse.cxx \
	$(NAVSRCDIR)/parse.h \
	$(NAVSRCDIR)/navsetup.cxx

bin_PROGRAMS = navcalc navsetup

navcalc_SOURCES = $(NAVHEADERS) $(COMMONSRC) $(CALCSRC)
navsetup_SOURCES = $(NAVHEADERS) $(COMMONSRC) $(SETUPSRC) $(TOKENIZERSRC)
