import numpy as np
import sympy as sym
import time as time
from math import exp,sqrt,sin,cos,pi
from scitools.std import *

from solvers import scalar, vectorized

'''
run a simple gaussian wave over the x axis. can be tested with both versions
'''
V = lambda x, y: 0
I = lambda x, y: 0
f = lambda x, y,t: 0
q = lambda x, y: 1

def I(x, y):
    return exp(-(x-10/2.0)**2/2.0 -(y-10/2.0)**2/2.0)


Lx = 10
Ly = 10

Nx = 30
Ny = 30

#vectorized(I, V, f, q, b=0, Lx=Lx, Nx=Nx, Ly=Ly, Ny=Ny, dt=min(float(Lx)/(Nx*2), float(Ly)/(Ny*2)) , T=10, plot=True)
scalar(I, V, f, q, b=0, Lx=Lx, Nx=Nx, Ly=Ly, Ny=Ny, dt=min(float(Lx)/(Nx*2), float(Ly)/(Ny*2)) , T=10, plot=True)
