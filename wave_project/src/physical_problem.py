import numpy as np
import sympy as sym
import time as time
from math import exp,sqrt,sin,cos,pi
from scitools.std import *

from solvers import scalar, vectorized

'''
This program simulates a wave through the sea with tree different shapes for hills
on the seafloor.
Usage: just run one of the functions defined below
'''

def ocean_waves(hill, version='vectorized'):
    #initial conditions
    V = lambda x, y: 0
    f = lambda x, y,t: 0
    def I(x, y):
        return 2*exp(-(x)**2/10.0)

    if hill == 'gaussian':
        def q(x,y):
            return 11 - 10*exp(-(x-30)**2/50.0 -(y-40)**2/50.0)
    if hill == 'cosine':
        def q(x,y):
            Mx = 30
            My = 40
            R = 15
            d = sqrt((x-Mx)**2+(y-My)**2)
            if d <=R:
                return 11 - 10*cos(pi*(x-Mx)/(2*R))*cos(pi*(y-My)/(2*R))
            else:
                return 11
    if hill == 'square':
        def q(x,y):
            if x <= 40 and x>=20 and y <=50 and y >= 30:
                return 1
            else:
                return 11
    
    #mesh parameters
    Lx = 80
    Ly = 80
    Nx = 100
    Ny = 100
    
    if version == 'scalar':
        scalar(I, V, f, q, b=0, Lx=Lx, Nx=Nx, Ly=Ly, Ny=Ny, dt=0.1 , T=40, plot=True)
    if version == 'vectorized':    
        vectorized(I, V, f, q, b=0, Lx=Lx, Nx=Nx, Ly=Ly, Ny=Ny, dt=0.1 , T=40, plot=True)


ocean_waves('gaussian')
