#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "lidDrivenCavity.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: Benny Malengier <bm@cage.ugent.be>
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  
 # ###################################################################
 ##

r"""

This example is an implementation of the well known lid driven cavity problem
at Re=1000, with comparison to the results in the literature. 
For this, the Navier-Stokes equations are solved on a colocated grid. 

This example is an extension of :mod:`examples.flow.stokesCavity` where the
viscous limit of this problem is considered. See that example for an 
explenation on the SIMPLE algorithm.

The equation that is solved is 

.. math::

   \vec{u} \cdot \nabla \vec{u} - \nabla \mu \cdot \nabla \vec{u} = -\nabla p
   
and the continuity equation,

.. math::
    
   \nabla \cdot \vec{u} = 0 
   
where :math:`\vec{u}` is the fluid velocity, :math:`p` is the pressure and :math:`\mu`
is the viscosity.  The domain in this example is a square cavity
of unit dimensions with a moving lid of unit speed.  This example
uses the SIMPLE algorithm with Rhie-Chow interpolation for colocated grids to solve
the pressure-momentum coupling, see also :mod:`examples.flow.stokesCavity`. 
A number of aspects of :term:`FiPy` need to be improved to have a first class flow
solver. These include, higher order spatial diffusion terms,
proper wall boundary conditions, improved mass flux evaluation and
extrapolation of cell values to the boundaries using gradients.

We compare the results with Ghia et al (1982), which is common in literature.
For an open source CFD solver that can solve the same problem, see Gerris_ , the values of Ghia et al where
taken from there 

.. _Gerris:    http://gfs.sourceforge.net/tests/tests/lid.html

To start, some parameters are declared.

>>> from fipy import *

>>> L = 1.0
>>> # number of cells, multiplicator of 3, not of 2, so 0.5 is a cell center!
>>> N = 63
>>> dL = L / N
>>> viscosity = 1e-3
>>> U = 1.
>>> Reynolds = U*L/viscosity
>>> #0.8 for pressure and 0.5 for velocity are typical relaxation values for SIMPLE
>>> pressureRelaxation = 0.8
>>> velocityRelaxation = 0.5
>>> if __name__ == '__main__':
...     sweeps = 500
... else:
...     sweeps = 15

Build the mesh.

.. index:: Grid2D
   
>>> mesh = Grid2D(nx=N, ny=N, dx=dL, dy=dL)

Declare the variables.

.. index:: CellVariable
   
>>> pressure = CellVariable(mesh=mesh, name='pressure')
>>> pressureCorrection = CellVariable(mesh=mesh)
>>> xVelocity = CellVariable(mesh=mesh, name='X velocity')
>>> yVelocity = CellVariable(mesh=mesh, name='Y velocity')

The velocity is required as a rank-1
:class:`~fipy.variables.faceVariable.FaceVariable` for calculating the mass
flux. This is required by the Rhie-Chow correction to avoid pressure/velocity
decoupling.

>>> velocity = FaceVariable(mesh=mesh, rank=1)
>>> velocity[0, mesh.getFacesTop().getValue()] = U

Build the Navier-Stokes equations in the cell centers. The previously calculated
velocity is taken as the coefficient of a convectionterm. As convectionterm we
choose a PowerLawConvectionTerm.

>>> xVelocityEq = PowerLawConvectionTerm(coeff=velocity) - \
...     DiffusionTerm(coeff=viscosity) + pressure.getGrad().dot([1.,0.])
>>> yVelocityEq = PowerLawConvectionTerm(coeff=velocity) - \
...     DiffusionTerm(coeff=viscosity) + pressure.getGrad().dot([0.,1.])
>>> #xVelocityEq = 2.* xVelocity * PowerLawConvectionTerm(coeff=(1.,0.)) \
... #                + yVelocity * PowerLawConvectionTerm(coeff=(0.,1.)) \
... #                + yVelocity.getGrad().dot([0.,1.])*xVelocity \
... #    - DiffusionTerm(coeff=viscosity) + pressure.getGrad().dot([1.,0.])
>>> #yVelocityEq =  2.* yVelocity * PowerLawConvectionTerm(coeff=(0.,1.)) \
... #                + xVelocity * PowerLawConvectionTerm(coeff=(1.,0.)) \
... #                + xVelocity.getGrad().dot([1.,0.])*yVelocity \
... #    - DiffusionTerm(coeff=viscosity) + pressure.getGrad().dot([0.,1.])

In this example the SIMPLE algorithm is used to couple the
pressure and momentum equations. We obtain the pressure correction equation, 

.. math::
    
   \nabla \frac{V_P}{a_P} \cdot \nabla p' = \nabla \cdot \vec{u}^{\ast}
   
In the discretized version of the above equation :math:`V_P / a_P` is
approximated at the face by :math:`A_f d_{AP} / (a_P)_f`. In :term:`FiPy` the
pressure correction equation can be written as, 

>>> ap = CellVariable(mesh=mesh, value=1.)
>>> coeff = 1./ ap.getArithmeticFaceValue() *mesh._getFaceAreas() * mesh._getCellDistances() 
>>> pressureCorrectionEq = DiffusionTerm(coeff=coeff) - velocity.getDivergence()

On a colocated grid as :term:`FiPy` uses, it is important to correctly define
:term:`velocity`, to avoid pressure oscillations. We apply the Rhie-Chow 
correction:

.. math::
    
  \vec{u}_f = \frac{1}{2}(\vec{u}^{\ast}_L + \vec{u}^{\ast}_R)) 
  + \frac{1}{2}\left(\frac{V}{a_P}\right)_{\mathrm{avg\ L,R}} (\nabla p^{\ast}_L+ \nabla p^{\ast}_R)
  - \left(\frac{V}{a_P}\right)_{\mathrm{avg\ L,R}} (\nabla p^{\ast}_f) 

where f is the face, and L and R the adjacent cells. We start by introducing needed
terms 

>>> from fipy.variables.faceGradVariable import _FaceGradVariable
>>> volume = CellVariable(mesh=mesh, value=mesh.getCellVolumes(), name='Volume')
>>> contrvolume=volume.getArithmeticFaceValue()

And set up the velocity with this formula in the SIMPLE loop. 
Now, set up the no-slip boundary conditions

.. index:: FixedValue
   
>>> bcs = (FixedValue(faces=mesh.getFacesLeft(), value=0.),
...        FixedValue(faces=mesh.getFacesRight(), value=0.),
...        FixedValue(faces=mesh.getFacesBottom(), value=0.),)
>>> bcsX = bcs + (FixedValue(faces=mesh.getFacesTop(), value=U),)
>>> bcsY = bcs + (FixedValue(faces=mesh.getFacesTop(), value=0.),)
>>> bcsPC = (FixedValue(faces=mesh.getFacesLeft() & (mesh.getFaceCenters()[1]<0.9*dL), value=0.),)

Set up the viewers,

.. index::
   :module: fipy.viewers
   
>>> if __name__ == '__main__':
...     viewer = Viewer(vars=(pressure, yVelocity, velocity, xVelocity),
...                xmin=0., xmax=1., ymin=0., ymax=1., colorbar=True)

Below, we iterate for a set number of sweeps. We use the :meth:`sweep`
method instead of :meth:`solve` because we require the residual for
output.  We also use the :meth:`cacheMatrix`, :meth:`getMatrix`,
:meth:`cacheRHSvector` and :meth:`getRHSvector` because both the matrix and
RHS vector are required by the SIMPLE algorithm. Additionally, the
:meth:`sweep` method is passed an ``underRelaxation`` factor to relax the
solution. This argument cannot be passed to :meth:`solve`.

>>> #solverpc = LinearCGSSolver(tolerance=1e-10)

.. index:: sweep, cacheMatrix, getMatrix, cacheRHSvector, getRHSvector
   
>>> for sweep in range(sweeps):
...
...     ## solve the Stokes equations to get starred values
...     xVelocityEq.cacheMatrix()
...     xres = xVelocityEq.sweep(var=xVelocity,
...                              boundaryConditions=bcsX,
...                              underRelaxation=velocityRelaxation)
...     yres = yVelocityEq.sweep(var=yVelocity,
...                              boundaryConditions=bcsY,
...                              underRelaxation=velocityRelaxation)
...     xmat = xVelocityEq.getMatrix()
...
...     ## update the ap coefficient from the matrix diagonal
...     ap[:] = xmat.takeDiagonal()
...
...     ## update the face velocities based on starred values with the 
...     ## Rhie-Chow correction. 
...     xvface = xVelocity.getArithmeticFaceValue()
...     yvface = yVelocity.getArithmeticFaceValue()
...     ## cell pressure gradient
...     presgrad = pressure.getGrad()
...     ## face pressure gradient
...     facepresgrad = _FaceGradVariable(pressure)
...
...     velocity[0] = xVelocity.getArithmeticFaceValue() \
...          + contrvolume / ap.getArithmeticFaceValue() * \
...            (presgrad[0].getArithmeticFaceValue()-facepresgrad[0])
...     velocity[1] = yVelocity.getArithmeticFaceValue() \
...          + contrvolume / ap.getArithmeticFaceValue() * \
...            (presgrad[1].getArithmeticFaceValue()-facepresgrad[1])
...     velocity[..., mesh.getExteriorFaces().getValue()] = 0.
...     velocity[0, mesh.getFacesTop().getValue()] = U
...
...     ## solve the pressure correction equation
...     pressureCorrectionEq.cacheRHSvector()
...     ## left bottom point must remain at pressure 0, so no correction
...     pres = pressureCorrectionEq.sweep(var=pressureCorrection, 
...                                       boundaryConditions=bcsPC,
...                                       #solver=solverpc
...                                       )
...     rhs = pressureCorrectionEq.getRHSvector()
...
...     ## update the pressure using the corrected value
...     pressure.setValue(pressure + pressureRelaxation * pressureCorrection)
...     #pressure.setValue(pressure + pressureRelaxation * (pressureCorrection -pressureCorrection.getGlobalValue()[0]))
...     ## update the velocity using the corrected pressure
...     xVelocity.setValue(xVelocity - pressureCorrection.getGrad()[0] / \
...                                                ap * mesh.getCellVolumes())
...     yVelocity.setValue(yVelocity - pressureCorrection.getGrad()[1] / \
...                                                ap * mesh.getCellVolumes())
...
...     if __name__ == '__main__':
...         if sweep%10 == 0:
...             print 'sweep:',sweep,', x residual:',xres, \
...                                  ', y residual',yres, \
...                                  ', p residual:',pres, \
...                                  ', continuity:',max(abs(rhs))
...
...             viewer.plot()

.. image:: lidcavityvx.*
   :width: 90%
   :align: center

.. image:: lidcavityvy.*
   :width: 90%
   :align: center

.. image:: lidcavity.*
   :width: 90%
   :align: center

Test values with the literature

>>> import numpy as np
>>> ## values of Ghia et al (1982)
>>> yexct =  np.array([-0.327052, -0.397406, -0.217948, -0.428914, -0.43629, -0.444335,
...     -0.046595, 0.001598, -0.4993, 0.118733, 0.235193, 0.352315, 0.45404, 0.461386,
...     0.469392, 0.476719, 0.5], float) + 0.5
>>> vxexact = [-0.383699, -0.297251, -0.27788, -0.222276, -0.201989, -0.181701,
...     -0.106804, -0.060949,  -0.000882, 0.057217, 0.186849, 0.333239, 0.466401,
...     0.511382, 0.574884, 0.659554, 0.999118]
>>> xexct = np.array([-0.500577, -0.43768, -0.429602, -0.421523, -0.406521,
...     -0.343624, -0.273803, -0.265724, -0.000289, 0.304962, 0.359781, 0.40652,
...     0.445182, 0.45326, 0.461339, 0.46884, 0.5],float) + 0.5
>>> vyexact = [0.00069404, 0.275621, 0.290847, 0.303994, 0.326826, 0.371038,
...     0.330015, 0.32307, 0.0252893, -0.318994, -0.427191, -0.515279, -0.392034,
...     -0.336623, -0.277749, -0.214023, -6.20706e-17]
>>> # calc values at y=0.5
>>> xcalc = mesh.getCellCenters()[0][N*(N//2):N*(N//2+1)]
>>> vycalc = yVelocity.getValue()[N*(N//2):N*(N//2+1)]
>>> # calc values at x=0.5
>>> ycalc = mesh.getCellCenters()[1][N//2::N]
>>> vxcalc = xVelocity.getValue()[N//2::N]
>>> print vxcalc[N//2], numerix.allclose(vxcalc[N//2], -0.0241933937349)
0
>>> print vycalc[N//2], numerix.allclose(vycalc[N//2], 0.00163975227638)
1
>>> ## profiles along y=0.5 of v_y and along x=0.5 of v_x
>>> if __name__ == '__main__':
...     import pylab
...     pylab.figure()
...     pylab.axis(xmin=0,xmax=1,ymin=np.min(vyexact)-0.02, ymax=np.max(vyexact)+0.02)
...     pylab.title('v_y velocity along y=%s compared with y=0.5 Ghia et al (1982)' % mesh.getCellCenters()[1][N*(N//2)])
...     pylab.plot(xcalc, vycalc, 'b-', xexct, vyexact, 'ro')
...     pylab.figure()
...     pylab.axis(xmin=0,xmax=1,ymin=np.min(vxexact)-0.02, ymax=np.max(vxexact)+0.02)
...     pylab.title('v_x velocity along x=%s compared with x=0.5 Ghia et al (1982)' % mesh.getCellCenters()[0][N//2])
...     pylab.plot(ycalc, vxcalc, 'b-', yexct, vxexact, 'ro')


.. image:: lidcavitylitcmpvx.*
   :width: 90%
   :align: center

.. image:: lidcavitylitcmpvy.*
   :width: 90%
   :align: center
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript(__name__))
    raw_input('finished')
