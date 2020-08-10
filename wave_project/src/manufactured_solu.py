import sympy as sym
import numpy as np
from solvers import *
import matplotlib.pyplot as plt
from matplotlib import cm

def source_term(w, mx, Lx, my, Ly, A, B, c, b):
    x, y, t =sym.symbols('x y t')
    kx = mx*np.pi / Lx
    ky = my*np.pi / Ly
    ue = (A*sym.cos(w*t) + B*sym.sin(w*t))*sym.exp(-c*t)*sym.cos(kx*x)*sym.cos(ky*y)
    q = 0.2
    utt = sym.diff(ue, t, t)
    ut = sym.diff(ue, t)
    uxq = sym.diff(q*sym.diff(ue,x),x)
    uyq = sym.diff(q*sym.diff(ue,y),y)
    
    st = sym.simplify(utt + b*ut - uxq - uyq)
    f = sym.lambdify((x,y,t), st)
    ute = sym.lambdify((x,y,t), ut)
    return f, ute

def test(w = 1, mx=1, Lx=10, my=1, Ly=10, A=1, B=1, b = 1, c=1):
    kx = mx*np.pi / Lx
    ky = my*np.pi / Ly
    def ue(x,y,t):
        return (A*np.cos(w*t) + B*np.sin(w*t))*np.exp(-c*t)*np.cos(kx*x)*np.cos(ky*y)
    I = lambda x, y : ue(x,y,0)
    q = lambda x, y : 0.2
    f, ut = source_term(w, mx, Lx, my, Ly, A, B, c, b)
    V = lambda x, y : ut(x,y,0)
    
    h = [0.25, 0.1, 0.07, 0.04]#, 0.1]
    rate = np.zeros(5)
    diff = np.zeros(5)

    
    # calculate convergence
    for i in range(4):
        dt = 0.1*h[i]
        Nx = int(round(float(Lx)/h[i]))
        Ny = int(round(float(Ly)/h[i]))
        diff[i] = scalar(I, V, f, q, b, Lx, Nx, Ly, Ny, dt, 1, ue, user_action='diff')
        if i > 0:
            rate[i] = np.log((diff[i]/diff[i-1]))/np.log((float(h[i])/h[i-1]))
        
    print Nx
    return rate, h, Nx

r,h,Nx = test()
print r, h
