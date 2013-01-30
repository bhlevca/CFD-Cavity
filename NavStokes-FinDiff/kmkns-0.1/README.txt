================================================================================
                                   kmkns-0.1
             A Navier-Stokes Solver for Single- and Two-Phase Flow
                       Created by Kim Motoyoshi Kalland
                                 September 2008
================================================================================

*** DESCRIPTION ***

kmkns is a 2D Navier-Stokes solver implemented as a Python package with Python 
modules and C++ extension modules. It uses the finite difference method on a 
uniform, rectangular grid. It handles single- and two-phase incompressible, 
Newtonian, laminar flow with obstacles.


*** INSTALLATION ***

This package requires the NumPy and Gnuplot modules and the SciTools package.

SciTools can be downloaded from: http://code.google.com/p/scitools/

Install kmkns by typing "python setup.py install" in the "kmkns-0.1" directory.


*** LICENSE ***

Released under the LGPL license. In addition, the SuperLU license applies to the
"poisson" module binary file and the SuperLU source files.


*** FILES ***

The following files have been written by me, Kim Motoyoshi Kalland:
kmkns/__init__.py
kmkns/csf.py
kmkns/Keyboard.py
kmkns/NavierStokes.py
poisson/poisson.cpp
poisson/sparse.h
scripts/advection_test2.py
scripts/advection_test.py
scripts/bubble_test2.py
scripts/bubble_test.py
scripts/channel_test2.py
scripts/channel_test.py
scripts/drop_test2.py
scripts/drop_test.py
scripts/ls_scripts.txt
scripts/poisson_test.py
scripts/pressure_test2.py
scripts/pressure_test.py
scripts/profile_test.py
scripts/rayleigh_taylor_test2.py
scripts/rayleigh_taylor_test.py
scripts/square_test2.py
scripts/square_test3.py
scripts/square_test.py
scripts/static_drop_test2.py
scripts/static_drop_test.py
scripts/step_test2.py
scripts/step_test.py
scripts/tutorial1.py
scripts/tutorial2.py
vof/main.cpp
vof/vof.cpp
vof/vof.h
README.txt
setup.py

The following files belong to the SuperLU library:
poisson/colamd.c
poisson/colamd.h
poisson/dasum.c
poisson/daxpy.c
poisson/dcolumn_bmod.c
poisson/dcolumn_dfs.c
poisson/dcopy.c
poisson/dcopy_to_ucol.c
poisson/dgscon.c
poisson/dgsequ.c
poisson/dgsrfs.c
poisson/dgssvx.c
poisson/dgstrf.c
poisson/dgstrs.c
poisson/dlacon.c
poisson/dlamch.c
poisson/dlangs.c
poisson/dlaqgs.c
poisson/dmemory.c
poisson/dmyblas2.c
poisson/dpanel_bmod.c
poisson/dpanel_dfs.c
poisson/dpivotgrowth.c
poisson/dpivotL.c
poisson/dpruneL.c
poisson/dsnode_bmod.c
poisson/dsnode_dfs.c
poisson/dsp_blas2.c
poisson/dsp_blas3.c
poisson/dtrsv.c
poisson/dutil.c
poisson/f2c.h
poisson/get_perm_c.c
poisson/heap_relax_snode.c
poisson/idamax.c
poisson/lsame.c
poisson/memory.c
poisson/mmd.c
poisson/relax_snode.c
poisson/slu_Cnames.h
poisson/slu_ddefs.h
poisson/slu_util.h
poisson/sp_coletree.c
poisson/sp_ienv.c
poisson/sp_preorder.c
poisson/superlu_timer.c
poisson/supermatrix.h
poisson/util.c
poisson/xerbla.c
SuperLU_License.txt

