#!/usr/bin/env python

import math
from numpy import *

def smooth(phi, mask, delta):
    """
    Smooth 'phi' using a 3x3 stencil. The original 'phi' is not
    modified.
    
    The returned ghost and obstacle cell values are undefined.
    """

    #stencil = array([[1,4,1],[4,16,4],[1,4,1]], float)
    stencil = array([[1,2,1],[2,4,2],[1,2,1]], float)
    
    if (stencil.shape[0] & 1) == 0 or (stencil.shape[1] & 1) == 0:
        raise ValueError, 'Stencil must have odd number of elements.'
    
    m = mask & 1
    
    mid = array(stencil.shape, int) / 2
    
    sum_phi = zeros(phi.shape, float)   # sum of weight times phi
    sum_w = sum_phi.copy()              # sum of weights
    
    for i in xrange(stencil.shape[0]):
        x_offs = i - mid[0]
        x1 = max(0, x_offs) + 1
        x2 = phi.shape[0] + min(0, x_offs) - 1
        for j in xrange(stencil.shape[1]):
            y_offs = j - mid[1]
            y1 = max(0, y_offs) + 1
            y2 = phi.shape[1] + min(0, y_offs) - 1
    
            sum_phi[x1-x_offs:x2-x_offs,y1-y_offs:y2-y_offs] \
                += stencil[i,j] * phi[x1:x2,y1:y2] * m[x1:x2,y1:y2]
            sum_w[x1-x_offs:x2-x_offs,y1-y_offs:y2-y_offs] \
                += stencil[i,j] * m[x1:x2,y1:y2]
            
    return sum_phi / (sum_w + (sum_w == 0.0)) # avoid division by zero

    
def gradient(phi, mask, delta):
    """
    Calculate the face centered gradient of phi which is cell centered.
    Neumann boundary conditions are applied. 
    """
    dx, dy = delta
    
    gx = zeros((phi.shape[0] - 1, phi.shape[1]), float)
    gy = zeros((phi.shape[0], phi.shape[1] - 1), float)
    
    gx[1:-1,1:-1] = (phi[2:-1,1:-1] - phi[1:-2,1:-1]) / dx
    gy[1:-1,1:-1] = (phi[1:-1,2:-1] - phi[1:-1,1:-2]) / dy
    
    gx[1:-1,1:-1] *= mask[1:-2,1:-1] & mask[2:-1,1:-1] & 1
    gy[1:-1,1:-1] *= mask[1:-1,1:-2] & mask[1:-1,2:-1] & 1
    
    return gx, gy
    
    
def curvature(gx, gy, delta):
    """
    Calculate the curvature from the given gradient.
    """
    dx, dy = delta
    div = zeros((gy.shape[0], gx.shape[1]), float)
    
    # Vertex centered gradients
    vgx = 0.5 * (gx[:,1:] + gx[:,:-1])
    vgy = 0.5 * (gy[1:,:] + gy[:-1,:])
    
    len_vg = (vgx**2 + vgy**2)**0.5
    
    div_g = (vgx[1:,1:] + vgx[1:,:-1] - vgx[:-1,1:] - vgx[:-1,:-1]) / (2.0 * dx) \
        + (vgy[1:,1:] - vgy[1:,:-1] + vgy[:-1,1:] - vgy[:-1,:-1]) / (2.0 * dy)
    
    cgx = 0.25 * (vgx[1:,1:] + vgx[1:,:-1] + vgx[:-1,1:] + vgx[:-1,:-1])
    cgy = 0.25 * (vgy[1:,1:] + vgy[1:,:-1] + vgy[:-1,1:] + vgy[:-1,:-1])
    
    len_cg = (cgx**2 + cgy**2)**0.5
    len_cg += (len_cg == 0.0)
    
    n_dot_grad = (len_vg[1:,1:] + len_vg[1:,:-1] - len_vg[:-1,1:] - len_vg[:-1,:-1]) \
        / (2.0 * dx) * (cgx / len_cg) \
        + (len_vg[1:,1:] - len_vg[1:,:-1] + len_vg[:-1,1:] - len_vg[:-1,:-1]) \
        / (2.0 * dy) * (cgy / len_cg)
    
    div[1:-1,1:-1] = (n_dot_grad - div_g) / len_cg
    
    return div

def surface_tension(phi, mask, delta, sigma):
    """
    Calculate the face centered surface tension force as:
    
    -sigma * kappa * grad(phi)
    
    where sigma is the surface tension coefficient.
    'phi' should be smooth.
    """
    
    phi_kappa = smooth(phi, mask, delta)
    gx, gy = gradient(phi_kappa, mask, delta)
    kappa = curvature(gx, gy, delta)

    gx, gy = gradient(phi, mask, delta)
    
    Fx = 0.5 * sigma * gx * (kappa[1:,:] + kappa[:-1,:])
    Fy = 0.5 * sigma * gy * (kappa[:,1:] + kappa[:,:-1])
    
    return Fx, Fy
    
    