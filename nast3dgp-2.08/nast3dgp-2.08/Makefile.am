## process this file with automake to produce Makefile.in
# This is the top level Makefile for NaSt3DGP.

SHELL = /bin/sh
@SET_MAKE@
ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src mpi doc tools
EXTRA_DIST = m4/acx_mpi.m4

# support for old targets:
setup:
	@cd src; $(MAKE) navsetup

serial: setup
	@cd src; $(MAKE) navcalc

parallel: setup
	@cd mpi; $(MAKE) all-am

.PHONY: setup serial parallel

nobase_dist_pkgdata_DATA = \
	examples/Cavity/cavity.nav \
	examples/RayleighBenard/rb.nav \
	examples/RayleighBenard/rb_glyc.nav \
	examples/RayleighBenard/rb3D.nav \
	examples/BackwardFacingStep/bwf_step.nav \
	tools/GridGen/GridGen.f \
	tools/vrml2nav/vrml2nav.cfg \
	tools/vrml2nav/examples/porsche_911.wrl \
	tools/vrml2nav/examples/submarine.wrl \
	tools/vrml2nav/examples/vrml2nav.cfg \
	tools/visualization/matlab/ReadNaSt3D.m \
	tools/visualization/matlab/StreamFunction.m \
	tools/visualization/matlab/cavity.m \
	tools/visualization/vtk/viewNaSt3DGP.tcl \
	tools/visualization/vtk/README

dist-hook:
	rm -rf `find $(distdir)/doc -name src`
