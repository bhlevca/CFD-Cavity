#!/usr/bin/env python

from numpy import *
import scitools.easyviz as viz
from kmkns.poisson import *
import math
import kmkns.vof
import kmkns.csf
import os

def _lerp(a, b, s):
    return a * (1 - s) + b * s


def _neumann_bc(phi, mask):
    """
    Extrapolate values from the fluid domain to ghost cells
    and obstacle cells.
    """
    w = (mask & 1) # weight = 1 for fluid cells, 0 for obstacle cells
    w_sum = zeros(w.shape, int) # sum of weights for each cell

    # set ghost cells to be treated as obstacle cells
    w[0,:] = 0
    w[-1,:] = 0
    w[:,0] = 0
    w[:,-1] = 0

    # first set obstacle values to zero, keep fluid values
    phi *= w

    # add right fluid cell to left obstacle cell
    weight = (1 - w[:-1,:]) * w[1:,:]
    phi[:-1,:] += weight * phi[1:,:]
    w_sum[:-1,:] += weight

    # add left fluid cell to right obstacle cell
    weight = (1 - w[1:,:]) * w[:-1,:]
    phi[1:,:] += weight * phi[:-1,:]
    w_sum[1:,:] += weight

    # add top fluid cell to bottom obstacle cell
    weight = (1 - w[:,:-1]) * w[:,1:]
    phi[:,:-1] += weight * phi[:,1:]
    w_sum[:,:-1] += weight

    # add bottom fluid cell to top obstacle cell
    weight = (1 - w[:,1:]) * w[:,:-1]
    phi[:,1:] += weight * phi[:,:-1]
    w_sum[:,1:] += weight

    # add top right fluid cell to bottom left obstacle cell
    # if it doesn't have a value already (sum of weights = 0).
    weight = (1 - w[:-1,:-1]) * (w_sum[:-1,:-1] == 0) * w[1:,1:]
    phi[:-1,:-1] += weight * phi[1:,1:]
    w_sum[:-1,:-1] += weight

    # add bottom left fluid cell to top right obstacle cell
    weight = (1 - w[1:,1:]) * (w_sum[1:,1:] == 0) * w[:-1,:-1]
    phi[1:,1:] += weight * phi[:-1,:-1]
    w_sum[1:,1:] += weight

    # add bottom right fluid cell to top left obstacle cell
    weight = (1 - w[:-1,1:]) * (w_sum[:-1,1:] == 0) * w[1:,:-1]
    phi[:-1,1:] += weight * phi[1:,:-1]
    w_sum[:-1,1:] += weight

    # add top left fluid cell to bottom right obstacle cell
    weight = (1 - w[1:,:-1]) * (w_sum[1:,:-1] == 0) * w[:-1,1:]
    phi[1:,:-1] += weight * phi[:-1,1:]
    w_sum[1:,:-1] += weight

    # divide by the sum of weights
    phi /= (w_sum + (w_sum == 0)) # avoid division by zero

    return phi


class Grid(object):
    def __init__(self, min_coord, max_coord, divs):
        """
        Create a uniform, rectangular grid.

        min_coord = lower, left corner of the grid
        max_coord = upper, right corner of the grid
        divs      = number of cells horizontally and vertically
        """
        super(Grid, self).__init__()
        self.min_coord = array(min_coord, float)
        self.max_coord = array(max_coord, float)
        self.divs = array(divs, int)


    def __repr__(self):
        return "Grid(" + ", ".join([repr(self.min_coord), repr(self.max_coord), repr(self.divs)]) + ")"


    def __getitem__(self, index):
        """
        Return node (x_i, y_j) where 'index' = (i,j).
        (x_0, y_0) is the lower left corner of the grid.
        (x_m, y_n) is the upper right corner where (m, n) == grid.shape.
        'i' and 'j' don't need to be integers and they don't need to
        be in range [0,m]x[0,n].
        """
        i = array(index, float)
        return self.min_coord + i/self.divs * (self.max_coord - self.min_coord)


    def get_shape(self):
        """
        Return the number of nodes (cell vertices) horizontally and 
        vertically.
        """
        return tuple(self.divs + 1)


    def get_delta(self):
        """
        Return the node spacing horizontally and vertically.
        """
        return (self.max_coord - self.min_coord)/self.divs


    def get_coordinate(self, i):
        """
        Return one of the coordinate components for all nodes.
        For instance, get_coordinate(1)[3,4] will return the
        y-coordinates for node (3,4), i.e. y_4.

        (x_0, y_0) is the lower left corner of the grid.
        (x_m, y_n) is the upper right corner where (m, n) == grid.shape.
        """
        x = zeros(self.divs + 1, float)
        lin = linspace(self.min_coord[i], self.max_coord[i], self.divs[i] + 1)
        s = (slice(None),)
        n = (len(self.divs) - 1 - i)
        for j in xrange(len(lin)):
            x[s*i + (j,) + s * n] = lin[j]
        return x

    shape = property(fget = get_shape)
    delta = property(fget = get_delta)
    x = property(fget = lambda self: self.get_coordinate(0))
    y = property(fget = lambda self: self.get_coordinate(1))
    z = property(fget = lambda self: self.get_coordinate(2))


class NavierStokes2D(object):
    EAST = 'east'
    WEST = 'west'
    NORTH = 'north'
    SOUTH = 'south'

    MASS_SCALE = 'scale'
    MASS_ADD = 'add'
    MASS_IGNORE = 'ignore'
    
    CURV_CSF = 'csf'
    CURV_DAC = 'dac'
    CURV_MDAC = 'mdac'

#========================#
# Initialization methods #
#========================#

    def __init__(self, grid, **kwargs):
        """
        Create a simple 2D Navier-Stokes solver.
        (Newtonian, incompressible, laminar, immiscible two-phase flow in a
        uniform grid)

        Possible arguments:
        mu or mu_liquid       -- Dynamic viscosity of main fluid.
        nu or nu_fluid        -- Kinematic viscosity of main fluid.
        rho or rho_liquid     -- Density of main fluid.
        mu_gas                -- Dynamic viscosity of second fluid.
        nu_gas                -- Kinematic viscosity of second fluid.
        rho_gas               -- Density of second fluid.
        surface_tension_coeff -- Surface tension coefficient.
        gravity               -- Acceleration due to gravity, must have two
                                 components.
        curv_method           -- Method used for curvature calculation:
                                 "csf"    = Continuum surface force method.
                                 "dac"    = Direction averaged curvature.
                                 "mdac"   = Modified DAC.
        mass_conservation     -- Method used for outflow correction:
                                 "scale"  = Scale outflow by constant factor.
                                 "add"    = Add constant to outflow.
                                 "ignore" = Leave the outflow as it is.
        property_smoothing    -- Smoothing of density and viscosity near 
                                 the interface:
                                 True     = On.
                                 False    = Off.

        Either mu or nu should be given, not both.
        """

        # Call superclass constructor. Probably not needed.
        #super(NavierStokes2D, self).__init__()

        # Initialize data structure.
        self._grid = grid       # stores some grid info
        self.t = 0.0
        self.p = zeros(array(grid.shape) + 1, float)
        self.u = zeros((self.p.shape[0] - 1, self.p.shape[1]), float)
        self.v = zeros((self.p.shape[0], self.p.shape[1] - 1), float)
        self.mask = ones(self.p.shape, int)
        self.c = zeros(self.p.shape, float) # vof-value
        self.dye = zeros(self.p.shape, float)
        self._bc = {self.EAST:{}, self.WEST:{}, self.NORTH:{}, self.SOUTH:{}}
        self._bc_finalized = False
        self._particles = []
        self.it = 0
        self._poisson_data = poisson_data()
        self._gas_cell_fraction = [] # the sum of (c > 0.5) divided by total area
        self._gas_fraction = [] # the sum of c divided by total area
        self._gas_centre = []
        self._liquid_centre = []
        self._record_time = []  # the time corresponding to the values stored in the lists above

        # Initialize parameters
        # Default parameters
        default_params = {'rho_liquid' : 1.0, 'rho_gas' : 0.0013,
                          'mu_liquid' : 1.0, 'mu_gas' : 0.016,
                          'surface_tension_coeff' : 0.0, 'gravity' : [0.0, 0.0],
                          'two_phase' : False, 'use_dye' : False,
                          'curv_method' : self.CURV_MDAC,
                          'mass_conservation' : self.MASS_ADD,
                          'property_smoothing' : False}
        for key in default_params:
            self.__dict__['_' + key] = default_params[key]

        # Merge optional parameter names
        if 'mu' in kwargs: kwargs['mu_liquid'] = kwargs['mu']
        if 'nu' in kwargs: kwargs['nu_fluid'] = kwargs['nu']
        if 'rho' in kwargs: kwargs['rho_liquid'] = kwargs['rho']

        if 'rho_liquid' in kwargs: self._rho_liquid = kwargs['rho_liquid']
        if 'rho_gas' in kwargs: self._rho_gas = kwargs['rho_gas']

        # Calculate dynamic viscosity
        if 'nu_fluid' in kwargs:
            kwargs['mu_liquid'] = kwargs['nu_fluid'] * self._rho_liquid
        if 'nu_gas' in kwargs:
            kwargs['mu_gas'] = kwargs['nu_gas'] * self._rho_gas

        # Store parameters
        for key in default_params:
            if key in kwargs:
                self.__dict__['_' + key] = kwargs[key]


    def set_ic(self, **kwargs):
        """
        Set initial condition for the domain. The keywords should be 'u', 'v',
        'p' or 'c'. The corresponding values must be either a constant or a
        function of 'x' and 'y' (space).
        """
        for ic in kwargs:
            # Skip invalid initial conditions. (Maybe this should raise an exception.)
            if not ic in ('u', 'v', 'p', 'c'):
                continue

            f = kwargs[ic]

            # Find the correct array to initialize.
            # If 'ic' is "u", 'a' is set to 'self.u', etc.
            a = self.__dict__[ic]

            # If 'f' is a function, call it on each node and store the result.
            if callable(f):
                # Calculate node offsets.
                # True is converted to integer one, False to integer zero.
                xOffset = 0.5 * (ic in ('p', 'v', 'c'))
                yOffset = 0.5 * (ic in ('p', 'u', 'c'))

                # iterate over all nodes.
                for i in xrange(a.shape[0]):
                    for j in xrange(a.shape[1]):
                        node = self._grid[i - xOffset, j - yOffset]
                        a[i,j] = f(node[0], node[1])
            else:
                a[:,:] = f # 'f' is a constant.


    def set_bc(self, boundaries, **kwargs):
        """
        Set boundary condition for the four boundaries east, west, north or
        south. 'boundaries' should be one of the strings 'east', 'west',
        'north' or 'south', or a combination with plus signs in between, eg.
        'west+north+south'. The keywords should be 'u', 'v', 'p', 'dudn',
        'dvdn' or 'c'. The corresponding values must be either a constant
        or a function of 'x', 'y' and 't' (space and time). Boundary
        conditions must be set before calling 'step()'.
        """
        for boundary in boundaries.split('+'):
            # If boundary is 'north', 'south', 'east' or 'west'...
            # (ignore otherwise)
            if boundary in self._bc:
                # Store boundary condition data.
                self._bc[boundary].update(kwargs)

                if 'p' in kwargs:
                    if boundary == self.NORTH:
                        self.mask[:, -1] = 3
                    elif boundary == self.SOUTH:
                        self.mask[:, 0] = 3
                    elif boundary == self.WEST:
                        self.mask[0, :] = 3
                    elif boundary == self.EAST:
                        self.mask[-1, :] = 3
                if 'v' in kwargs:
                    if boundary in (self.NORTH, self.SOUTH):
                        self._bc[boundary]['dpdn'] = 0
                if 'u' in kwargs:
                    if boundary in (self.EAST, self.WEST):
                        self._bc[boundary]['dpdn'] = 0
                # u, v, dpdn -> mask = 1, which is the default


    def _remove_singularity(self):
        """
        If Neumann boundary condition is used along the entire boundary,
        the Poisson equation becomes singular. This function fixes one or
        more of the fluid pressure nodes if necessary such that the system
        has one unique solution.
        """
        remove_singularity(self.mask) # defined in the poisson module


    def reset_time(self):
        self.t = 0.0


    def add_obstacle(self, min_coord, max_coord):
        """
        Add impenetrable cells to the mask.
        Cells whose centre lies within the given rectangle
        are marked as impenetrable.

        min_coord = lower left corner
        max_coord = upper right corner

        The coordinates must be inside the domain or on the boundary.
        """
        scale = self._grid.divs/(self._grid.max_coord - self._grid.min_coord)
        # Calculate lower and upper bound in number of cells.
        # Add one for the ghost cells.
        low = scale * (array(min_coord, float) - self._grid.min_coord) + 1
        high = scale * (array(max_coord, float) - self._grid.min_coord) + 1
        self.mask[round(low[0]):round(high[0]),round(low[1]):round(high[1])] = 0


    def set_obstacles(self, f):
        """
        Initialize mask with obstacles.
        'f' is called for each cell with the coordinate of the cell's 
        centre as parameter. If 'f' evaluates to 'True', the
        cell is an obstacle cell. Otherwise, the cell is a fluid cell.

        f = function of x and y, returns True for obstacle cells
        """
        mask = self.mask
        for i in xrange(1, mask.shape[0] - 1):
            for j in xrange(1, mask.shape[1] - 1):
                node = self._grid[i - 0.5, j - 0.5]
                if f(node[0], node[1]):
                    mask[i,j] = 0
                else: 
                    mask[i,j] = 1
                    

    def add_particle(self, pos):
        """
        Add a particle to the fluid. The particle will follow the fluid flow.
        However, in the current implementation, particles tend to fall out of
        the domain in the corners. Adding particles to the fluid will slow
        down the simulation considerably.
        """
        self._particles.append(array(pos, float))


    def add_tracer(self, min_coord, max_coord):
        """
        Add a block of dye to the fluid. The dye will follow the fluid flow
        without diffusion.

        min_coord = lower left corner
        max_coord = upper right corner

        The coordinates must be inside the domain or on the boundary.
        """
        self._use_dye = True
        low = self._grid.min_coord
        high = self._grid.max_coord
        delta = self._grid.delta
        x = linspace(low[0] - delta[0], high[0] + delta[0], self._grid.divs[0] + 3)
        y = linspace(low[1] - delta[1], high[1] + delta[1], self._grid.divs[1] + 3)
        self.dye = maximum(self.dye,
            kmkns.vof.rectangle_fraction(min_coord[0], min_coord[1],
            max_coord[0], max_coord[1], x, y))


    def add_gas(self, min_coord, max_coord):
        """
        Mark a block of fluid as second fluid.

        min_coord = lower left corner
        max_coord = upper right corner
        """
        self._two_phase = True
        low = self._grid.min_coord
        high = self._grid.max_coord
        delta = self._grid.delta
        x = linspace(low[0] - delta[0], high[0] + delta[0], self._grid.divs[0] + 3)
        y = linspace(low[1] - delta[1], high[1] + delta[1], self._grid.divs[1] + 3)
        self.c = maximum(self.c,
            kmkns.vof.rectangle_fraction(min_coord[0], min_coord[1],
            max_coord[0], max_coord[1], x, y))


    def add_bubble(self, centre, radius, scale_x=1.0, scale_y=1.0):
        """
        Mark a circular area of fluid as second fluid.
        """
        self._two_phase = True
        low = self._grid.min_coord
        high = self._grid.max_coord
        delta = self._grid.delta
        x = linspace(low[0] - delta[0], high[0] + delta[0], self._grid.divs[0] + 3) / scale_x
        y = linspace(low[1] - delta[1], high[1] + delta[1], self._grid.divs[1] + 3) / scale_y
        self.c = maximum(self.c,
            kmkns.vof.circle_fraction(centre[0] / scale_x,
            centre[1] / scale_y, radius, x, y))


    def add_liquid(self, min_coord, max_coord):
        """
        Mark a block of fluid as main fluid.

        min_coord = lower left corner
        max_coord = upper right corner
        """
        self._two_phase = True
        low = self._grid.min_coord
        high = self._grid.max_coord
        delta = self._grid.delta
        x = linspace(low[0] - delta[0], high[0] + delta[0], self._grid.divs[0] + 3)
        y = linspace(low[1] - delta[1], high[1] + delta[1], self._grid.divs[1] + 3)
        self.c = minimum(self.c,
            1.0 - kmkns.vof.rectangle_fraction(min_coord[0], min_coord[1],
            max_coord[0], max_coord[1], x, y))


    def add_drop(self, centre, radius, scale_x=1.0, scale_y=1.0):
        """
        Mark a circular area of fluid as main fluid.
        """
        self._two_phase = True
        low = self._grid.min_coord
        high = self._grid.max_coord
        delta = self._grid.delta
        x = linspace(low[0] - delta[0], high[0] + delta[0], self._grid.divs[0] + 3) / scale_x
        y = linspace(low[1] - delta[1], high[1] + delta[1], self._grid.divs[1] + 3) / scale_y
        self.c = minimum(self.c,
            1.0 - kmkns.vof.circle_fraction(centre[0] / scale_x,
            centre[1] / scale_y, radius, x, y))

    def set_interface(self, orientation, F):
        """
        Create an interface between liquid and gas.
        
        orientation = orientation of the interface:
                      0 = liquid to the south
                      1 = liquid to the east
                      2 = liquid to the north
                      3 = liquid to the west
        F           = integral of the liquid depth function.
                      If orientation is 0 or 2, F must be a function of x.
                      If orientation is 1 or 2, F must be a function of y.
        """
        self._two_phase = True
        low, high = self._grid.min_coord, self._grid.max_coord
        delta, c = self._grid.delta, self.c
        x = linspace(low[0] - delta[0], high[0] + delta[0], self._grid.divs[0] + 3)
        y = linspace(low[1] - delta[1], high[1] + delta[1], self._grid.divs[1] + 3)
        
        if orientation in (0, 2):
            for i in xrange(len(x) - 1):
                h = (F(x[i+1]) - F(x[i])) / (delta[0] * delta[1]) + 1.0
                n = floor(h) # number of liquid cells
                r = h - n # remaining liquid
                if orientation == 0:
                    c[i,:n] = 0.0
                    c[i,n:] = 1.0
                    c[i,n] -= r
                else:
                    c[i,-n:] = 0.0
                    c[i,:-n] = 1.0
                    c[i,-n-1] -= r
        else:
            for i in xrange(len(y) - 1):
                h = (F(y[i+1]) - F(y[i])) / (delta[0] * delta[1]) + 1.0
                n = floor(h) # number of liquid cells
                r = h - floor(h) # remaining liquid
                if orientation == 3:
                    c[:n,i] = 0.0
                    c[n:,i] = 1.0
                    c[n,i] -= r
                else:
                    c[-n:,i] = 0.0
                    c[:-n,i] = 1.0
                    c[-n-1,i] -= r
        

#====================#
# Simulation methods #
#====================#

    def step(self, dt, beta=1.0, gamma=0.0):
        """
        Advance the simulation to the next time c.

        dt    = time stepping value (delta time)
        beta  = factor in the projection method
                See "H. P. Langtangen, K.-A. Mardal and R. Winther:
                Numerical Methods for Incompressible Viscous Flow, 2002"
        gamma = Upwind differencing factor
                See "Numerical Simulation in Fluid Dynamics: A Practical 
                Introduction" (1997) by Griebel, Dornsheifer and Neunhoeffer.
        """

        # Retrieve node spacing and flow variables
        delta, mask, c = self._grid.delta, self.mask, self.c
        u, v, p = self.u, self.v, self.p

        if not self._bc_finalized:
            self._remove_singularity()
            self._bc_finalized = True

        # Impose boundary conditions.
        self._update_velocity_bc()
        self._update_pressure_bc()
        self._update_level_bc()

        # Update particles.
        for i in xrange(len(self._particles)):
            self._particles[i] += dt * self._interpolate_velocity(self._particles[i])

        # If dye was added, update the dye.
        if self._use_dye:
            kmkns.vof.advect(self.dye, u, v, mask, delta, dt, self.it)

        # Calculate density and viscosity
        if self._two_phase:
            # kmkns.vof.advect() modifies 'c' in place, thus 'self.c' is
            # updated.
            kmkns.vof.advect(c, u, v, mask, delta, dt, self.it)
            _neumann_bc(c, mask)

            # smooth does not modify 'c' in place, thus
            # 'self.c' is not modified.
            if self._property_smoothing:
                c = kmkns.csf.smooth(c, mask, delta)
                _neumann_bc(c, mask)

        rho = self._rho_liquid * (1.0 - c) + self._rho_gas * c
        mu = self._mu_liquid * (1.0 - c) + self._mu_gas * c

        rho_u = 0.5 * (rho[1:,:] + rho[:-1,:])
        rho_v = 0.5 * (rho[:,1:] + rho[:,:-1])

        # Calculate tentative (predicted) velocity.
        (du, dv) = self._calc_velocity_change(dt, rho, rho_u, rho_v, mu, beta, gamma)

        u += du
        v += dv

        # Update time.
        self.t += dt
        self.it += 1

        # Impose BC to tentative velocity.
        self._update_tentative_velocity_bc()

        f = self._poisson_rhs(dt)

        # Calculate boundary conditions for 'phi' and store the values in
        # the ghost cells.
        phi = zeros(p.shape, float)
        self._update_phi_bc(phi, beta)

        # Calculate pressure correction. (phi = p^(n+1) - beta * p^(n))
        # 'mask' defines where to use Dirichlet and Neumann boundary conditions.
        p *= beta
        if self._two_phase:
            p += poisson(phi, f, delta, mask, rho, self._poisson_data, 1) # re-factorise
        else:
            p += poisson(phi, f, delta, mask, rho, self._poisson_data, 3) # re-use LU

        # Clear pressure inside obstacles.
        p[1:-1,1:-1] *= (mask[1:-1,1:-1] & 1)

        # Correct velocity. Zero velcity change on obstacle walls.
        u[1:-1,1:-1] -= dt * (phi[2:-1,1:-1] - phi[1:-2,1:-1]) / \
            (delta[0] * rho_u[1:-1,1:-1]) * (mask[2:-1,1:-1] & mask[1:-2,1:-1] & 1)
        v[1:-1,1:-1] -= dt * (phi[1:-1,2:-1] - phi[1:-1,1:-2]) / \
            (delta[1] * rho_v[1:-1,1:-1]) * (mask[1:-1,2:-1] & mask[1:-1,1:-2] & 1)

        # Correct velocity on boundary for Dirichlet pressure boundary condition
        u[0,1:-1] -= dt * (phi[1,1:-1] - phi[0,1:-1]) / \
            (delta[0] * rho_u[0,1:-1]) * (mask[1,1:-1] & 1)
        u[-1,1:-1] -= dt * (phi[-1,1:-1] - phi[-2,1:-1]) / \
            (delta[0] * rho_u[-1,1:-1]) * (mask[-2,1:-1] & 1)
        v[1:-1,0] -= dt * (phi[1:-1,1] - phi[1:-1,0]) / \
            (delta[1] * rho_v[1:-1,0]) * (mask[1:-1,1] & 1)
        v[1:-1,-1] -= dt * (phi[1:-1,-1] - phi[1:-1,-2]) / \
            (delta[1] * rho_v[1:-1,-1]) * (mask[1:-1,-2] & 1)

        if self._two_phase:
            self._update_records(rho)


    def _poisson_rhs(self, dt):
        """
        Calculate the right hand side of the Poisson equation
        'div(rho * grad(phi)) = -f'. Make sure the equation is
        solvable by using the method specified by the
        '_mass_conservation' attribute.
        """

        # f[1:-1,1:-1].sum() must be zero unless Dirichlet
        # pressure boundary condition is used.

        u, v, p, delta, mask = self.u, self.v, self.p, self._grid.delta, self.mask

        net_outflow = (u[-1,1:-1] - u[0,1:-1]).sum() * delta[1] + \
            (v[1:-1,-1] - v[1:-1,0]).sum() * delta[0]
        outflow_length = 0.0
        outflow = 0.0
        dirichlet_used = ('p' in self._bc[self.NORTH]) or \
            ('p' in self._bc[self.SOUTH]) or \
            ('p' in self._bc[self.WEST]) or \
            ('p' in self._bc[self.EAST])

        if 'dvdn' in self._bc[self.NORTH]:
            outflow_length += mask[1:-1,-2].sum() * delta[0]
            outflow += v[1:-1,-1].sum() * delta[0]
        if 'dvdn' in self._bc[self.SOUTH]:
            outflow_length += mask[1:-1,1].sum() * delta[0]
            outflow -= v[1:-1,0].sum() * delta[0]
        if 'dudn' in self._bc[self.WEST]:
            outflow_length += mask[1,1:-1].sum() * delta[1]
            outflow -= u[0,1:-1].sum() * delta[1]
        if 'dudn' in self._bc[self.EAST]:
            outflow_length += mask[-2,1:-1].sum() * delta[1]
            outflow += u[-1,1:-1].sum() * delta[1]

        if (not dirichlet_used) and outflow_length > 0.0 and \
            self._mass_conservation in (self.MASS_ADD, self.MASS_SCALE):

            if outflow == 0.0 or self._mass_conservation == self.MASS_ADD:
                flow_corr = net_outflow / outflow_length
                if 'dvdn' in self._bc[self.NORTH]:
                    v[1:-1,-1] -= mask[1:-1,-2] * flow_corr
                if 'dvdn' in self._bc[self.SOUTH]:
                    v[1:-1,0] += mask[1:-1,1] * flow_corr
                if 'dudn' in self._bc[self.WEST]:
                    u[0,1:-1] += mask[1,1:-1] * flow_corr
                if 'dudn' in self._bc[self.EAST]:
                    u[-1,1:-1] -= mask[-2,1:-1] * flow_corr
            else:
                flow_corr = 1.0 - net_outflow / outflow
                if 'dvdn' in self._bc[self.NORTH]:
                    v[1:-1,-1] *= flow_corr
                if 'dvdn' in self._bc[self.SOUTH]:
                    v[1:-1,0] *= flow_corr
                if 'dudn' in self._bc[self.WEST]:
                    u[0,1:-1] *= flow_corr
                if 'dudn' in self._bc[self.EAST]:
                    u[-1,1:-1] *= flow_corr

        # Calculate right hand side of the Poisson equation "-phi''=f"
        f = zeros(p.shape, float)
        f[1:-1,1:-1] = ((u[1:,1:-1]-u[:-1,1:-1]) / delta[0] + \
                        (v[1:-1,1:]-v[1:-1,:-1]) / delta[1])

        # third option for mass conservation
        if (not dirichlet_used) and (outflow_length == 0.0 or \
            not self._mass_conservation in (self.MASS_ADD, self.MASS_SCALE)):
            f[1:-1,1:-1] -= f[1:-1,1:-1].sum() / f[1:-1,1:-1].size

        return (-1.0 / dt) * f


    def _calc_velocity_change(self, dt, rho, rho_u, rho_v, mu, beta=1.0, gamma=0.0):
        """
        Calculate the tentative minus the previous velocity (u^*-u^n, v^*-u^n)
        If 'beta' is 1, this becomes the acceleration times 'dt'.
        """
        delta, mask = self._grid.delta, self.mask
        u, v, p = self.u, self.v, self.p

        du = zeros(self.u.shape, float)
        dv = zeros(self.v.shape, float)

        mu_off = 0.25 * (mu[1:,1:] + mu[1:,:-1] + mu[:-1,1:] + mu[:-1,:-1])


        # Viscosity
        # d/dy(mu*(du/dy+dv/dx))/rho
        du[:,1:-1] += (mu_off[:,1:] * ((u[:,2:] - u[:,1:-1]) / (delta[1]**2) + \
                      (v[1:,1:] - v[:-1,1:]) / (delta[0] * delta[1])) - \
                      mu_off[:,:-1] * ((u[:,1:-1] - u[:,:-2]) / (delta[1]**2) + \
                      (v[1:,:-1] - v[:-1,:-1]) / (delta[0]*delta[1]))) / \
                      rho_u[:,1:-1]
        # d/dx(2*mu*du/dx)/rho
        du[1:-1,1:-1] += (mu[2:-1,1:-1] * (u[2:,1:-1] - u[1:-1,1:-1]) - \
                          mu[1:-2,1:-1] * (u[1:-1,1:-1] - u[:-2,1:-1])) / \
                         (0.5 * (delta[0]**2) * rho_u[1:-1,1:-1])
        # boundary (assuming dirichlet bc for pressure)
        # (if another bc is used, these values will be overwritten later in the algorithm)
        du[0,1:-1] -= (mu[1,1:-1] * (v[1,1:] - v[1,:-1]) - \
                       mu[0,1:-1] * (v[0,1:] - v[0,:-1])) / \
                      (0.5 * delta[0] * delta[1] * rho_u[0,1:-1])
        du[-1,1:-1] -= (mu[-1,1:-1] * (v[-1,1:] - v[-1,:-1]) - \
                        mu[-2,1:-1] * (v[-2,1:] - v[-2,:-1])) / \
                       (0.5 * delta[0] * delta[1] * rho_u[-1,1:-1])

        # d/dx(mu*(du/dy+dv/dx))/rho
        dv[1:-1,:] += (mu_off[1:,:] * ((v[2:,:] - v[1:-1,:])/(delta[0]**2) + \
                      (u[1:,1:] - u[1:,:-1]) / (delta[0]*delta[1])) - \
                      mu_off[:-1,:] * ((v[1:-1,:] - v[:-2,:]) / (delta[0]**2) + \
                      (u[:-1,1:] - u[:-1,:-1])/(delta[0]*delta[1]))) / \
                      rho_v[1:-1,:]
        # d/dy(2*mu*dv/dy)/rho
        dv[1:-1,1:-1] += (mu[1:-1,2:-1] * (v[1:-1,2:] - v[1:-1,1:-1]) - \
                          mu[1:-1,1:-2] * (v[1:-1,1:-1] - v[1:-1,:-2])) / \
                         (0.5 * (delta[1]**2) * rho_v[1:-1,1:-1])
        # boundary
        dv[1:-1,0] -= (mu[1:-1,1] * (u[1:,1] - u[:-1,1]) - \
                       mu[1:-1,0] * (u[1:,0] - u[:-1,0])) / \
                      (0.5 * delta[0] * delta[1] * rho_v[1:-1,0])
        dv[1:-1,-1] -= (mu[1:-1,-1] * (u[1:,-1] - u[:-1,-1]) - \
                        mu[1:-1,-2] * (u[1:,-2] - u[:-1,-2])) / \
                       (0.5 * delta[0] * delta[1] * rho_v[1:-1,-1])


        # Convection
        # d(uu)/dx
        du[1:-1,1:-1] -= 0.25 * ((u[1:-1,1:-1] + u[2:,1:-1])**2 - \
                                 (u[1:-1,1:-1] + u[:-2,1:-1])**2) / delta[0]
        # boundary 2u*du/dx => -2u*dv/dy
        du[0,1:-1] += u[0,1:-1] * ((v[0,1:] + v[1,1:]) - \
                                   (v[0,:-1] + v[1,:-1])) / delta[1]
        du[-1,1:-1] += u[-1,1:-1] * ((v[-2,1:] + v[-1,1:]) - \
                                     (v[-2,:-1] + v[-1,:-1])) / delta[1]
        # d(uv)/dy
        du[:,1:-1] -= 0.25 * ((u[:,1:-1] + u[:,2:]) * (v[:-1,1:]+v[1:,1:]) - \
            (u[:,1:-1] + u[:,:-2]) * (v[:-1,:-1] + v[1:,:-1])) / delta[1]

        # d(vv)/dy
        dv[1:-1,1:-1] -= 0.25 * ((v[1:-1,1:-1] + v[1:-1,2:])**2 - \
                                 (v[1:-1,1:-1] + v[1:-1,:-2])**2) / delta[1]
        # boundary
        dv[1:-1,0] += v[1:-1,0] * ((u[1:,0] + u[1:,1]) - \
                                   (u[:-1,0] + u[:-1,1])) / delta[0]
        dv[1:-1,-1] += v[1:-1,-1] * ((u[1:,-2] + u[1:,-1]) - \
                                     (u[:-1,-2] + u[:-1,-1])) / delta[0]
        # d(vu)/dx
        dv[1:-1,:] -= 0.25 * ((v[1:-1,:] + v[2:,:]) * (u[1:,:-1] + u[1:,1:]) - \
            (v[1:-1,:] + v[:-2,:]) * (u[:-1,:-1] + u[:-1,1:])) / delta[0]

        # Upwind difference for convection (ref. griebel p. 29)
        # d(uu)/dx
        du[1:-1,1:-1] -= gamma * 0.25 * (abs(u[1:-1,1:-1] + u[2:,1:-1]) * \
            (u[1:-1,1:-1] - u[2:,1:-1]) - abs(u[:-2,1:-1]+u[1:-1,1:-1]) * \
            (u[:-2,1:-1] - u[1:-1,1:-1])) / delta[0]
        # d(uv)/dy
        du[1:-1,1:-1] -= gamma * 0.25 * ((u[1:-1,1:-1] - u[1:-1,2:]) * \
            abs(v[1:-2,1:] + v[2:-1,1:]) - (u[1:-1,:-2] - u[1:-1,1:-1]) * \
            abs(v[1:-2,:-1] + v[2:-1,:-1])) / delta[1]

        # d(vv)/dy
        dv[1:-1,1:-1] -= gamma * 0.25 * (abs(v[1:-1,1:-1] + v[1:-1,2:]) * \
            (v[1:-1,1:-1] - v[1:-1,2:]) - abs(v[1:-1,1:-1] + v[1:-1,:-2]) * \
            (v[1:-1,:-2] - v[1:-1,1:-1])) / delta[1]
        # d(vu)/dx
        dv[1:-1,1:-1] -= gamma * 0.25 * ((v[1:-1,1:-1] - v[2:,1:-1]) * \
            abs(u[1:,1:-2] + u[1:,2:-1]) - (v[:-2,1:-1] - v[1:-1,1:-1]) * \
            abs(u[:-1,1:-2] + u[:-1,2:-1])) / delta[0]


        # gravity
        du += self._gravity[0]
        dv += self._gravity[1]


        # surface tension
        if self._two_phase and self._surface_tension_coeff != 0.0:
            F_x, F_y = self._surface_tension()
            F_x /= rho_u[1:-1,1:-1]
            F_y /= rho_v[1:-1,1:-1]
            du[1:-1,1:-1] += F_x
            dv[1:-1,1:-1] += F_y


        # pressure forces
        # dp/dx
        du[:,1:-1] -= (beta / delta[0]) * (p[1:,1:-1] - p[:-1,1:-1]) / rho_u[:,1:-1]
        # dp/dy
        dv[1:-1,:] -= (beta / delta[1]) * (p[1:-1,1:] - p[1:-1,:-1]) / rho_v[1:-1,:]

        return (du * dt, dv * dt)


    def _surface_tension(self):
        c = self.c
        delta = self._grid.delta
        mask = self.mask
        u, v = self.u, self.v

        if self._curv_method == self.CURV_CSF:
            F_x, F_y = kmkns.csf.surface_tension(self.c, mask, delta,
                self._surface_tension_coeff)
            return F_x[1:-1,1:-1], F_y[1:-1,1:-1]

        #c = _neumann_bc(c.copy(), mask)

        kappa = kmkns.vof.curvature(c, mask, delta,
            self._curv_method == self.CURV_MDAC)[1:-1,1:-1]

        # Calculate weights
        w = c[1:-1,1:-1] * (1.0 - c[1:-1,1:-1]) + 1e-16

        grad_x = (c[1:,1:-1] - c[:-1,1:-1]) / delta[0]
        grad_y = (c[1:-1,1:] - c[1:-1,:-1]) / delta[1]

        # Impose Neumann BC for obstacles
        grad_x[1:-1,:] *= (mask[2:-1,1:-1] & mask[1:-2,1:-1] & 1)
        grad_y[:,1:-1] *= (mask[1:-1,2:-1] & mask[1:-1,1:-2] & 1)

        F_x = self._surface_tension_coeff * (kappa[1:,:] * w[1:,:] + \
            kappa[:-1,:] * w[:-1,:]) * grad_x[1:-1,:] / (w[1:,:] + w[:-1,:])
        F_y = self._surface_tension_coeff * (kappa[:,1:] * w[:,1:] + \
            kappa[:,:-1] * w[:,:-1]) * grad_y[:,1:-1] / (w[:,1:] + w[:,:-1])

        return F_x, F_y


    def _update_tentative_velocity_bc(self):
        """
        Update the velocity for prescribed velocity boundary condition.
        """
        u, v = self.u, self.v

        # Interior boundaries
        # Apply no-slip boundary conditions to obstacles.
        # Setup masks that are 0 where velocities need to be updated,
        # and 1 where they stay unmodified.
        # Note that (mask & 1) has 1 in the ghost cells.
        u_mask = (self.mask[:-1,:] | self.mask[1:,:]) & 1
        v_mask = (self.mask[:,:-1] | self.mask[:,1:]) & 1

        # zero velocity inside and on the boundary of obstacles
        u[:,:] *= (self.mask[:-1,:] & self.mask[1:,:] & 1)
        # negate velocities inside obstacles
        u[:,1:-2] -= (1 - u_mask[:,1:-2]) * u[:,2:-1]
        u[:,2:-1]  -= (1 - u_mask[:,2:-1]) * u[:,1:-2]

        # zero velocity inside and on the boundary of obstacles
        v[:,:] *= (self.mask[:,:-1] & self.mask[:,1:] & 1)
        # negate velocities inside obstacles
        v[1:-2,:] -= (1 - v_mask[1:-2,:]) * v[2:-1,:]
        v[2:-1,:]  -= (1 - v_mask[2:-1,:]) * v[1:-2,:]

        # north boundary
        bc = self._bc[self.NORTH]
        if 'v' in bc:
            f = bc['v']
            if callable(f):
                for i in xrange(v.shape[0]):
                    node = self._grid[i-0.5, v.shape[1]-1]
                    v[i,-1] = f(node[0], node[1], self.t) * (self.mask[i,-2] & 1)
            else:
                v[:,-1] = f

        # south boundary
        bc = self._bc[self.SOUTH]
        if 'v' in bc:
            f = bc['v']
            if callable(f):
                for i in xrange(v.shape[0]):
                    node = self._grid[i-0.5, 0]
                    v[i,0] = f(node[0], node[1], self.t) * (self.mask[i,1] & 1)
            else:
                v[:,0] = f

        # east boundary
        bc = self._bc[self.EAST]
        if 'u' in bc:
            f = bc['u']
            if callable(f):
                for i in xrange(u.shape[1]):
                    node = self._grid[u.shape[0]-1, i-0.5]
                    u[-1,i] = f(node[0], node[1], self.t) * (self.mask[-2,i] & 1)
            else:
                u[-1,:] = f

        # west boundary
        bc = self._bc[self.WEST]
        if 'u' in bc:
            f = bc['u']
            if callable(f):
                for i in xrange(u.shape[1]):
                    node = self._grid[0, i-0.5]
                    u[0,i] = f(node[0], node[1], self.t) * (self.mask[1,i] & 1)
            else:
                u[0,:] = f
                

    def _update_velocity_bc(self):
        u, v, c = self.u, self.v, self.c

        # prescribed velocity boundary condition
        self._update_tentative_velocity_bc()

        # teh rest
        # north boundary
        bc = self._bc[self.NORTH]
        if 'dvdn' in bc:
            if bc['dvdn'] != 0:
                raise ValueError, 'dv/dn must be zero.'
            v[:,-1] = v[:,-2] # extrapolate
        # 'v' in bc was handled in _update_tentative_velocity_bc()
        
        if 'dudn' in bc:
            if bc['dudn'] != 0:
                raise ValueError, 'du/dn must be zero.'
            u[:,-1] = u[:,-2] # extrapolate
        elif 'u' in bc:
            f = bc['u']
            if callable(f):
                for i in xrange(u.shape[0]):
                    node = self._grid[i, u.shape[1]-2]
                    u[i,-1] = 2 * f(node[0], node[1], self.t) - u[i,-2]
            else:
                u[:,-1] = 2 * f - u[:,-2]

        # south boundary
        bc = self._bc[self.SOUTH]
        if 'dvdn' in bc:
            if bc['dvdn'] != 0:
                raise ValueError, 'dv/dn must be zero.'
            v[:,0] = v[:,1]
        
        if 'dudn' in bc:
            if bc['dudn'] != 0:
                raise ValueError, 'du/dn must be zero.'
            u[:,0] = u[:,1]
        elif 'u' in bc:
            f = bc['u']
            if callable(f):
                for i in xrange(u.shape[0]):
                    node = self._grid[i, 0]
                    u[i,0] = 2*f(node[0], node[1], self.t) - u[i,1]
            else:
                u[:,0] = 2*f - u[:,1]

        # east boundary
        bc = self._bc[self.EAST]
        if 'dudn' in bc:
            if bc['dudn'] != 0:
                raise ValueError, 'du/dn must be zero.'
            u[-1,:] = u[-2,:]
        
        if 'dvdn' in bc:
            if bc['dvdn'] != 0:
                raise ValueError, 'dv/dn must be zero.'
            v[-1,:] = v[-2,:]
        elif 'v' in bc:
            f = bc['v']
            if callable(f):
                for i in xrange(v.shape[1]):
                    node = self._grid[v.shape[0]-2, i]
                    v[-1,i] = 2*f(node[0], node[1], self.t) - v[-2,i]
            else:
                v[-1,:] = 2*f - v[-2,:]

        # west boundary
        bc = self._bc[self.WEST]
        if 'dudn' in bc:
            if bc['dudn'] != 0:
                raise ValueError, 'du/dn must be zero.'
            u[0,:] = u[1,:]
        
        if 'dvdn' in bc:
            if bc['dvdn'] != 0:
                raise ValueError, 'dv/dn must be zero.'
            v[0,:] = v[1,:]
        elif 'v' in bc:
            f = bc['v']
            if callable(f):
                for i in xrange(v.shape[1]):
                    node = self._grid[0, i]
                    v[0,i] = 2*f(node[0], node[1], self.t) - v[1,i]
            else:
                v[0,:] = 2*f - v[1,:]


    def _update_pressure_bc(self):
        p = self.p

        # north boundary
        bc = self._bc[self.NORTH]
        if 'dpdn' in bc:
            if bc['dpdn'] != 0:
                raise ValueError, 'dp/dn must be zero.'
            p[:,-1] = p[:,-2]
        elif 'p' in bc:
            f = bc['p']
            if callable(f):
                for i in xrange(p.shape[0]):
                    node = self._grid[i-0.5, p.shape[1]-2]
                    p[i,-1] = 2 * f(node[0], node[1], self.t) - p[i,-2]
            else:
                p[:,-1] = 2 * f - p[:,-2]

        # south boundary
        bc = self._bc[self.SOUTH]
        if 'dpdn' in bc:
            if bc['dpdn'] != 0:
                raise ValueError, 'dp/dn must be zero.'
            p[:,0] = p[:,1]
        elif 'p' in bc:
            f = bc['p']
            if callable(f):
                for i in xrange(p.shape[0]):
                    node = self._grid[i-0.5, 0]
                    p[i,0] = 2 * f(node[0], node[1], self.t) - p[i,1]
            else:
                p[:,0] = 2 * f - p[:,1]

        # east boundary
        bc = self._bc[self.EAST]
        if 'dpdn' in bc:
            if bc['dpdn'] != 0:
                raise ValueError, 'dp/dn must be zero.'
            p[-1,:] = p[-2,:]
        elif 'p' in bc:
            f = bc['p']
            if callable(f):
                for i in xrange(p.shape[1]):
                    node = self._grid[p.shape[0]-2, i-0.5]
                    p[-1,i] = 2 * f(node[0], node[1], self.t) - p[-2,i]
            else:
                p[-1,:] = 2 * f - p[-2,:]

        # east boundary
        bc = self._bc[self.WEST]
        if 'dpdn' in bc:
            if bc['dpdn'] != 0:
                raise ValueError, 'dp/dn must be zero.'
            p[0,:] = p[1,:]
        elif 'p' in bc:
            f = bc['p']
            if callable(f):
                for i in xrange(p.shape[1]):
                    node = self._grid[0, i-0.5]
                    p[0,i] = 2 * f(node[0], node[1], self.t) - p[1,i]
            else:
                p[0,:] = 2 * f - p[1,:]


    def _update_phi_bc(self, phi, beta):
        p = self.p

        # north boundary
        bc = self._bc[self.NORTH]
        if 'p' in bc:
            f = bc['p']
            if callable(f):
                for i in xrange(p.shape[0]):
                    node = self._grid[i-0.5, p.shape[1]-2]
                    phi[i,-1] = 2 * f(node[0], node[1], self.t) - beta * (p[i,-2] + p[i,-1])
            else:
                phi[:,-1] = 2 * f - beta * (p[:,-2] + p[:,-1])

        # south boundary
        bc = self._bc[self.SOUTH]
        if 'p' in bc:
            f = bc['p']
            if callable(f):
                for i in xrange(p.shape[0]):
                    node = self._grid[i-0.5, 0]
                    phi[i,0] = 2 * f(node[0], node[1], self.t) - beta * (p[i,1] + p[i,0])
            else:
                phi[:,0] = 2 * f - beta * (p[:,1] + p[:,0])

        # east boundary
        bc = self._bc[self.EAST]
        if 'p' in bc:
            f = bc['p']
            if callable(f):
                for i in xrange(p.shape[1]):
                    node = self._grid[p.shape[0]-2, i-0.5]
                    phi[-1,i] = 2 * f(node[0], node[1], self.t) - beta * (p[-2,i] + p[-1,i])
            else:
                phi[-1,:] = 2 * f - beta * (p[-2,:] + p[-1,:])

        # east boundary
        bc = self._bc[self.WEST]
        if 'p' in bc:
            f = bc['p']
            if callable(f):
                for i in xrange(p.shape[1]):
                    node = self._grid[0, i-0.5]
                    phi[0,i] = 2 * f(node[0], node[1], self.t) - beta * (p[1,i] + p[0,i])
            else:
                phi[0,:] = 2 * f - beta * (p[1,:] + p[0,:])


    def _update_level_bc(self):
        c = self.c

        # north boundary
        bc = self._bc[self.NORTH]
        if 'c' in bc:
            f = bc['c']
            if callable(f):
                for i in xrange(c.shape[0]):
                    node = self._grid[i - 0.5, c.shape[1] - 1.5]
                    c[i,-1] = f(node[0], node[1], self.t)
            else:
                c[:,-1] = f

        # south boundary
        bc = self._bc[self.SOUTH]
        if 'c' in bc:
            f = bc['c']
            if callable(f):
                for i in xrange(c.shape[0]):
                    node = self._grid[i - 0.5, -0.5]
                    c[i,0] = f(node[0], node[1], self.t)
            else:
                c[:,0] = f

        # east boundary
        bc = self._bc[self.EAST]
        if 'c' in bc:
            f = bc['c']
            if callable(f):
                for i in xrange(c.shape[1]):
                    node = self._grid[c.shape[0] - 1.5, i - 0.5]
                    c[-1,i] = f(node[0], node[1], self.t)
            else:
                c[-1,:] = f

        # east boundary
        bc = self._bc[self.WEST]
        if 'c' in bc:
            f = bc['c']
            if callable(f):
                for i in xrange(c.shape[1]):
                    node = self._grid[-0.5, i - 0.5]
                    c[0,i] = f(node[0], node[1], self.t)
            else:
                c[0,:] = f


    def _update_records(self, rho):
        """
        Updates the _gas_fraction, _record_time and _gas_centre lists.
        (Obstacle cells erroneously counts as fluid cells.)
        """
        if not self._two_phase:
            return
            
        # TODO: Obstacle cells should not count.
            
        c = self.c    
        self._record_time.append(self.t)
        s = c[1:-1,1:-1].size
        gas_sum = c[1:-1,1:-1].sum()
        liquid_sum = s - gas_sum
        
        self._gas_cell_fraction.append((c[1:-1,1:-1] > 0.5).sum() / float(s))
        self._gas_fraction.append(gas_sum / s)

        low, high = self._grid.min_coord, self._grid.max_coord
        delta = self._grid.delta
        x = linspace(low[0] + 0.5 * delta[0], high[0] - 0.5 * delta[0], self._grid.divs[0])
        y = linspace(low[1] + 0.5 * delta[1], high[1] - 0.5 * delta[1], self._grid.divs[1])

        x_mc = (c[1:-1,1:-1].transpose() * x).sum() / gas_sum
        y_mc = (c[1:-1,1:-1] * y).sum() / gas_sum
        self._gas_centre.append([x_mc, y_mc])

        _1_minus_c = 1 - c
        x_mc = (_1_minus_c[1:-1,1:-1].transpose() * x).sum() / liquid_sum
        y_mc = (_1_minus_c[1:-1,1:-1] * y).sum() / liquid_sum
        self._liquid_centre.append([x_mc, y_mc])
        

#=============================================#
# Methods for retrieving and visualizing data #
#=============================================#

    def _interpolate_u(self, pos):
        """
        Return the u-velocity at the given position in the domain by using
        bilinear interpolation.
        """
        pos = (pos - self._grid.min_coord)/(self._grid.max_coord - self._grid.min_coord)
        pos = array([pos[0] * self._grid.divs[0], pos[1] * self._grid.divs[1] + 0.5])
        low = array(pos, int)
        if low[0] < 0: low[0] = 0
        if low[0] > self.u.shape[0]-2: low[0] = self.u.shape[0]-2
        if low[1] < 0: low[1] = 0
        if low[1] > self.u.shape[1]-2: low[1] = self.u.shape[1]-2
        frac = pos - low
        left = _lerp(self.u[low[0], low[1]], self.u[low[0], low[1]+1], frac[1])
        right = _lerp(self.u[low[0]+1, low[1]], self.u[low[0]+1, low[1]+1], frac[1])
        return _lerp(left, right, frac[0])


    def _interpolate_v(self, pos):
        """
        Return the v-velocity at the given position in the domain by using
        bilinear interpolation.
        """
        pos = (pos - self._grid.min_coord)/(self._grid.max_coord - self._grid.min_coord)
        pos = array([pos[0] * self._grid.divs[0] + 0.5, pos[1] * self._grid.divs[1]])
        low = array(pos, int)
        if low[0] < 0: low[0] = 0
        if low[0] > self.v.shape[0]-2: low[0] = self.v.shape[0]-2
        if low[1] < 0: low[1] = 0
        if low[1] > self.v.shape[1]-2: low[1] = self.v.shape[1]-2
        frac = pos - low
        left = _lerp(self.v[low[0], low[1]], self.v[low[0], low[1]+1], frac[1])
        right = _lerp(self.v[low[0]+1, low[1]], self.v[low[0]+1, low[1]+1], frac[1])
        return _lerp(left, right, frac[0])


    def _interpolate_velocity(self, pos):
        return array([self._interpolate_u(pos), self._interpolate_v(pos)])


    def avg_velocity(self, middle=True):
        """
        Calculate the velocity in the middle of the cells using linear
        interpolation if 'middle' is true. If 'middle' is false, the
        average velocity is calculated on grid cell corners.
        """
        if middle:
            u_avg = 0.5 * (self.u[:-1,1:-1] + self.u[1:,1:-1])
            v_avg = 0.5 * (self.v[1:-1,:-1] + self.v[1:-1,1:])
        else:
            u_avg = 0.5 * (self.u[:,1:] + self.u[:,:-1])
            v_avg = 0.5 * (self.v[1:,:] + self.v[:-1,:])

        return (u_avg, v_avg)


    def plot_gas_fraction(self, fig=1):
        """
        Plot the area of the second fluid relative to
        total fluid.
        """
        viz.figure(fig)
        viz.plot(self._record_time, self._gas_fraction, 'b-',
            self._record_time, self._gas_cell_fraction, 'r-')


    def plot_liquid_fraction(self, fig=1):
        """
        Plot the area of the first fluid (c <= 1/2) relative to
        total fluid.
        """
        viz.figure(fig)
        viz.plot(self._record_time, [1.0 - x for x in self._gas_fraction], 'b-',
            self._record_time, [1.0 - x for x in self._gas_cell_fraction], 'r-')


    def _mass_centre(self):
        """
        Calculate the mass centre from self._gas_centre and 
        self._liquid_centre.
        """
        cm = [[self._rho_gas * G[0] + self._rho_liquid * L[0],
            self._rho_gas * G[1] + self._rho_liquid * L[1]] 
            for G, L in zip(self._gas_centre, self._liquid_centre)]
        return cm
        

    def plot_mass_centre(self, fig=1):
        """
        Plot the path of the mass centre.
        """
        cm = self._mass_centre()
        viz.figure(fig)
        x = [c[0] for c in cm]
        y = [c[1] for c in cm]
        viz.plot(x, y)


    def plot_gas_centre(self, fig=1):
        """
        Plot the path of the gas centre.
        """
        viz.figure(fig)
        x = [c[0] for c in self._gas_centre]
        y = [c[1] for c in self._gas_centre]
        viz.plot(x, y)


    def plot_liquid_centre(self, fig=1):
        """
        Plot the path of the liquid centre.
        """
        viz.figure(fig)
        x = [c[0] for c in self._liquid_centre]
        y = [c[1] for c in self._liquid_centre]
        viz.plot(x, y)


    def plot_mass_centre_velocity(self, fig=1):
        """
        Plot the absolute value of the velocity of the mass centre
        with respect to time.
        """
        viz.figure(fig)
        mc = array(self._mass_centre(), float)
        t = array(self._record_time, float)
        v = (mc[1:] - mc[:-1])
        v[:,0] /= (t[1:] - t[:-1])
        v[:,1] /= (t[1:] - t[:-1])
        if len(t) > 1:
            viz.plot(0.5*(t[1:] + t[:-1]), (v[:,0]**2 + v[:,1]**2)**0.5)


    def plot_gas_centre_velocity(self, fig=1):
        """
        Plot the absolute value of the velocity of the gas centre
        with respect to time.
        """
        viz.figure(fig)
        mc = array(self._gas_centre, float)
        t = array(self._record_time, float)
        v = (mc[1:] - mc[:-1])
        v[:,0] /= (t[1:] - t[:-1])
        v[:,1] /= (t[1:] - t[:-1])
        if len(t) > 1:
            viz.plot(0.5*(t[1:] + t[:-1]), (v[:,0]**2 + v[:,1]**2)**0.5)


    def plot_liquid_centre_velocity(self, fig=1):
        """
        Plot the absolute value of the velocity of the liquid centre
        with respect to time.
        """
        viz.figure(fig)
        mc = array(self._liquid_centre, float)
        t = array(self._record_time, float)
        v = (mc[1:] - mc[:-1])
        v[:,0] /= (t[1:] - t[:-1])
        v[:,1] /= (t[1:] - t[:-1])
        if len(t) > 1:
            viz.plot(0.5*(t[1:] + t[:-1]), (v[:,0]**2 + v[:,1]**2)**0.5)


    def plot_surface_tension(self, fig=1):
        """
        Plot the curvature of the isolines of the c set function.
        """
        viz.figure(fig)

        delta = self._grid.delta
        shape = (self.p.shape[0] - 2, self.p.shape[1] - 2)
        F_x = zeros(shape, float)
        F_y = zeros(shape, float)
        F_x[:-1,:], F_y[:,:-1] = self._surface_tension()
        F_x[1:,:] += F_x[:-1,:].copy() # copy() is absolutely necessary
        F_y[:,1:] += F_y[:,:-1].copy()
        #F_x *= 0.5
        #F_y *= 0.5

        # scale the velocity to fit the plot
        sqrMax = (F_x[:,:]**2 + F_y[:,:]**2).max()
        if(sqrMax == 0.0): sqrMax = 1.0
        scale = 2.0*sqrt((delta[0]**2+delta[1]**2)/sqrMax)
        F_x *= scale
        F_y *= scale

        x = linspace(self._grid.min_coord[0] + 0.5*delta[0], self._grid.max_coord[0] - 0.5*delta[0], self._grid.divs[0])
        y = linspace(self._grid.min_coord[1] + 0.5*delta[1], self._grid.max_coord[1] - 0.5*delta[1], self._grid.divs[1])
        # the automatic scaling is no good. turn it off by passing zero.
        viz.quiver(x, y, F_x, F_y, 0)


    def plot_tracer(self, fig=1):
        """
        Plot the tracer in the fluid.
        """
        viz.figure(fig)
        delta = self._grid.delta
        x = linspace(self._grid.min_coord[0] + 0.5*delta[0], self._grid.max_coord[0] - 0.5*delta[0], self._grid.divs[0])
        y = linspace(self._grid.min_coord[1] + 0.5*delta[1], self._grid.max_coord[1] - 0.5*delta[1], self._grid.divs[1])
        viz.pcolor(x, y, self.dye[1:-1,1:-1], shading='flat')


    def plot_particles(self, fig=1):
        """
        Plot the particles as dots.
        """
        if len(self._particles) == 0:
            return
        viz.figure(fig)
        x = [self._particles[i][0] for i in xrange(len(self._particles))]
        y = [self._particles[i][1] for i in xrange(len(self._particles))]
        viz.plot(x, y, 'ko', axis=[self._grid.min_coord[0], self._grid.max_coord[0],
            self._grid.min_coord[1], self._grid.max_coord[1]])


    def plot_velocity(self, fig=1):
        """
        Plot the velocity as a vector field.
        """
        viz.figure(fig)

        delta = self._grid.delta
        (u_avg, v_avg) = self.avg_velocity()

        # scale the velocity to fit the plot
        sqrMaxU = (u_avg[:,:]**2 + v_avg[:,:]**2).max()
        if(sqrMaxU == 0.0): sqrMaxU = 1.0
        scale = 2.0*sqrt((delta[0]**2+delta[1]**2)/sqrMaxU)
        u_avg *= scale
        v_avg *= scale

        x = linspace(self._grid.min_coord[0] + 0.5*delta[0], self._grid.max_coord[0] - 0.5*delta[0], self._grid.divs[0])
        y = linspace(self._grid.min_coord[1] + 0.5*delta[1], self._grid.max_coord[1] - 0.5*delta[1], self._grid.divs[1])
        # the automatic scaling is no good. turn it off by passing zero.
        viz.quiver(x, y, u_avg, v_avg, 0)


    def plot_velocity_component(self, component, contours=None, fig=1):
        """
        Plot one of the velocity components.
        """
        viz.figure(fig)

        delta = self._grid.delta
        u_avg = self.avg_velocity()[component]

        x = linspace(self._grid.min_coord[0] + 0.5*delta[0], self._grid.max_coord[0] - 0.5*delta[0], self._grid.divs[0])
        y = linspace(self._grid.min_coord[1] + 0.5*delta[1], self._grid.max_coord[1] - 0.5*delta[1], self._grid.divs[1])
        if contours == None:
            viz.pcolor(x, y, u_avg, shading='flat')
        else:
            viz.contour(x, y, u_avg, contours, clabels='on')


    def plot_velocity_profile(self, component, slices, fig=1):
        """
        Plot the profile of one of the velocity components.
        'slices' is the position of the velocity profile, or a list of such
        positions. The positions must be in ascending order.
        """
        viz.figure(fig)

        if not hasattr(slices, '__iter__'):
            slices = [slices]

        velocity = (self.u, self.v)[component]
        amp = max(velocity.max(), 0.0) - min(velocity.min(), 0.0)

        if len(slices) > 1:
            s = array(slices)
            dist = (s[1:]-s[:-1]).min()
        else:
            dist = amp

        scale = amp/dist

        delta = self._grid.delta
        n = linspace(self._grid.min_coord[component^1] - 0.5*delta[component^1],
            self._grid.max_coord[component^1] + 0.5*delta[component^1],
            self._grid.divs[component^1] + 2)
        n[0] = self._grid.min_coord[component^1]
        n[-1] = self._grid.max_coord[component^1]

        data = []
        legend = []
        for t in slices:
            s = self._grid.divs[component]*(t-self._grid.min_coord[component])/(self._grid.max_coord[component]-self._grid.min_coord[component])
            s_int = min(int(floor(s)), self._grid.divs[component] - 1)
            s_frac = s - s_int
            if component == 0:
                v_lerp = velocity[s_int,:]*(1.0-s_frac) + velocity[s_int+1,:]*s_frac
                v_lerp[0] = 0.5*(v_lerp[0] + v_lerp[1])
                v_lerp[-1] = 0.5*(v_lerp[-1] + v_lerp[-2])
                data.append(v_lerp + t*scale)
                data.append(n)
                data.append(ones(self._grid.divs[component^1] + 2) * t*scale)
                data.append(n)
                data.append('k-')
                legend.append('x = ' + str(t))
                legend.append('')
            if component == 1:
                v_lerp = velocity[:,s_int]*(1.0-s_frac) + velocity[:,s_int+1]*s_frac
                v_lerp[0] = 0.5*(v_lerp[0] + v_lerp[1])
                v_lerp[-1] = 0.5*(v_lerp[-1] + v_lerp[-2])
                data.append(n)
                data.append(v_lerp + t*scale)
                data.append(n)
                data.append(ones(self._grid.divs[component^1] + 2) * t*scale)
                data.append('k-')
                legend.append('y = ' + str(t))
                legend.append('')
        viz.plot(*data, **{'legend':legend})


    def plot_level(self, contours=None, fig=1):
        """
        Plot the c set field. The nodes are located at grid cell centres.
        """
        viz.figure(fig)

        delta = self._grid.delta

        x = linspace(self._grid.min_coord[0] + 0.5*delta[0], self._grid.max_coord[0] - 0.5*delta[0], self._grid.divs[0])
        y = linspace(self._grid.min_coord[1] + 0.5*delta[1], self._grid.max_coord[1] - 0.5*delta[1], self._grid.divs[1])
        if contours == None:
            viz.pcolor(x, y, self.c[1:-1,1:-1], shading='flat')
        else:
            viz.contour(x, y, self.c[1:-1,1:-1], contours, clabels='on')


    def plot_pressure(self, contours=None, fig=1):
        """
        Plot the pressure field. The nodes are located at grid cell centres.
        """
        viz.figure(fig)

        delta = self._grid.delta

        x = linspace(self._grid.min_coord[0] + 0.5*delta[0], self._grid.max_coord[0] - 0.5*delta[0], self._grid.divs[0])
        y = linspace(self._grid.min_coord[1] + 0.5*delta[1], self._grid.max_coord[1] - 0.5*delta[1], self._grid.divs[1])
        if contours == None:
            viz.pcolor(x, y, self.p[1:-1,1:-1], shading='flat')
        else:
            viz.contour(x, y, self.p[1:-1,1:-1], contours, clabels='on')


    def plot_curvature(self, contours=None, fig=1):
        viz.figure(fig)

        delta = self._grid.delta

        x = linspace(self._grid.min_coord[0] + 0.5*delta[0], self._grid.max_coord[0] - 0.5*delta[0], self._grid.divs[0])
        y = linspace(self._grid.min_coord[1] + 0.5*delta[1], self._grid.max_coord[1] - 0.5*delta[1], self._grid.divs[1])

        if self._curv_method == self.CURV_CSF:
            phi = kmkns.csf.smooth(self.c, self.mask, delta)
            gx, gy = kmkns.csf.gradient(phi, self.mask, delta)
            f = kmkns.csf.curvature(gx, gy, delta)
        else:    
            f = kmkns.vof.curvature(self.c, self.mask, self._grid.delta,
                self._curv_method == self.CURV_MDAC)

        if contours == None:
            viz.pcolor(x, y, f[1:-1,1:-1], shading='flat')
        else:
            viz.contour(x, y, f[1:-1,1:-1], contours, clabels='on')


    def divergence(self):
        """
        Return the divergence field. The nodes coincide with pressure nodes.
        """
        delta = self._grid.delta
        f = zeros(self.p.shape, float)
        u = self.u
        v = self.v
        f[1:-1,1:-1] = (u[1:,1:-1]-u[:-1,1:-1])/delta[0] + (v[1:-1,1:]-v[1:-1,:-1])/delta[1]
        return f


    def plot_divergence(self, contours=None, fig=1):
        """
        Plot the divergence field.
        """
        viz.figure(fig)

        delta = self._grid.delta
        x = linspace(self._grid.min_coord[0] + 0.5*delta[0], self._grid.max_coord[0] - 0.5*delta[0], self._grid.divs[0])
        y = linspace(self._grid.min_coord[1] + 0.5*delta[1], self._grid.max_coord[1] - 0.5*delta[1], self._grid.divs[1])
        f = self.divergence()
        if contours == None:
            viz.pcolor(x, y, f[1:-1,1:-1], shading='flat')
        else:
            viz.contour(x, y, f[1:-1,1:-1], contours, clabels='on')


    def vorticity(self):
        """
        Return the vorticity or curl field. The nodes are located
        at the grid cell corners.
        """
        delta = self._grid.delta
        u = self.u
        v = self.v
        f = (v[1:,:]-v[:-1,:])/delta[0] - (u[:,1:]-u[:,:-1])/delta[1]
        return f


    def plot_vorticity(self, contours=None, fig=1):
        """
        Plot the vorticity.
        """
        viz.figure(fig)

        delta = self._grid.delta
        x = linspace(self._grid.min_coord[0], self._grid.max_coord[0], self._grid.divs[0] + 1)
        y = linspace(self._grid.min_coord[1], self._grid.max_coord[1], self._grid.divs[1] + 1)
        f = self.vorticity()
        if contours == None:
            viz.pcolor(x, y, f, shading='flat')
        else:
            viz.contour(x, y, f, contours, clabels='on')


    def stream_function(self):
        """
        Return the stream function field. The nodes are located
        at the grid cell corners.
        """
        u = self.u
        v = self.v
        phi = zeros(self._grid.shape, float)
        mask = (self.mask[:-1,1:-1] | self.mask[1:,1:-1]) & 1
        phi[1:,0] = -(v[1:-1,0] * self._grid.delta[0]).cumsum(0)
        phi[:,1:] = mask * u[:,1:-1] * self._grid.delta[1]
        phi = phi.cumsum(1)
        return phi


    def plot_stream(self, contours=None, fig=1):
        """
        Plot the stream function.
        """
        viz.figure(fig)

        delta = self._grid.delta
        x = linspace(self._grid.min_coord[0], self._grid.max_coord[0], self._grid.divs[0] + 1)
        y = linspace(self._grid.min_coord[1], self._grid.max_coord[1], self._grid.divs[1] + 1)
        f = self.stream_function()
        if contours == None:
            viz.pcolor(x, y, f, shading='flat')
        else:
            viz.contour(x, y, f, contours, clabels='on')


    def find_suitable_dt(self, safety_factor=0.95):
        """
        Calculate a time stepping value that should result in a stable
        simulation. This should be called once every time step.
        """
        delta = self._grid.delta
        u = abs(self.u) + 0.000001 * (self.u == 0.0) # avoid division by zero
        v = abs(self.v) + 0.000001 * (self.v == 0.0)

        if self._two_phase:
            nus = array([self._mu_liquid / self._rho_liquid,
                self._mu_gas / self._rho_gas], float)

            dt1 = 0.5 * min(delta[0] / u.max(), delta[1] / v.max())
            dt2 = 0.5 / (nus.max() * (delta**-2).sum())
            dt3 = 2.0 * nus.min() / max(u.max(), v.max())**2
            if self._surface_tension_coeff > 0.0:
                dt4 = ((self._rho_liquid + self._rho_gas) * delta.min()**3 / \
                    (4.0 * math.pi * self._surface_tension_coeff))**0.5
                return min(dt1,dt2,dt3,dt4) * safety_factor
            else:
                return min(dt1,dt2,dt3) * safety_factor
        else:
            nu = self._mu_liquid / self._rho_liquid

            dt1 = min(delta[0] / u.max(), delta[1] / v.max())
            # if tracer fluid is used, make sure to fulfill the time-step 
            # restriction for the volume-of-fluid advection algorithm.
            if self._use_dye:
                dt1 *= 0.5
            dt2 = 0.5 / (nu * (delta**-2).sum())
            dt3 = 2.0 * nu / max(u.max(), v.max())**2
            return min(dt1, dt2, dt3) * safety_factor


    def load(self, filename):
        """
        Load the flow variables from the given file. The size of the arrays
        in the file and in 'self' must match. 
        
        Return true if file exsists, false if not found.
        Exceptions are thrown for other errors.
        """
        
        if not os.path.exists(filename):
            return False
        infile = open(filename, "r")
        data = infile.read()
        infile.close()
        # Avoid using 'eval' so that only files with the correct format is loaded.
        # Remove commas, put spaces around '[' and ']' and split by whitespace
        split_data = " ".join(" ] ".join(" [ ".join(data.split("[")).split("]")).split(",")).split()
        def build(tokens, i):
            l = []
            while i < len(tokens):
                if tokens[i] == '[':
                    sub, i = build(tokens, i+1)
                    l.append(sub)
                elif tokens[i] == ']':
                    return (l, i)
                else:
                    l.append(float(tokens[i]))
                i += 1
            return (l, i)
        l, i = build(split_data, 0)
        if len(l) < 11:
            raise IndexError, "Incorrect file format: '%s'" % filename

        self.p = array(l[0], float)
        self.u = array(l[1], float)
        self.v = array(l[2], float)
        self.c = array(l[3], float)
        self.dye = array(l[4], float)
        self._particles = [array(e,float) for e in l[5]]
        self._gas_centre = l[6]
        self._liquid_centre = l[7]
        self._record_time = l[8]
        self._gas_cell_fraction = l[9]
        self._gas_fraction = l[10]
        self.t = l[11]
        self.it = int(l[12])

        return True


    def save(self, filename):
        """
        Save the flow variables to the given file.
        """
        outfile = open(filename, "w")
        outfile.write(repr(self.p.tolist()) + '\n' + \
            repr(self.u.tolist()) + '\n' + \
            repr(self.v.tolist()) + '\n' + \
            repr(self.c.tolist()) + '\n' + \
            repr(self.dye.tolist()) + '\n' + \
            repr([e.tolist() for e in self._particles]) + '\n' + \
            repr(self._gas_centre) + '\n' + \
            repr(self._liquid_centre) + '\n' + \
            repr(self._record_time) + '\n' + \
            repr(self._gas_cell_fraction) + '\n' + \
            repr(self._gas_fraction) + '\n' + \
            repr(self.t) + '\n' + \
            repr(self.it))
        outfile.close()
