from solvers import *
import numpy as np

def standing_undamped(A, mx, my):
    
    # constants    
    Lx = 1
    kx = mx*np.pi/Lx
    Ly = 1
    ky = my*np.pi/Ly
    b = 0   
    T = 0.1
    w = np.sqrt(2)*np.pi

    # functions
    I = lambda x, y: A*np.cos(kx*x)*np.cos(kx*y)
    q = lambda x, y: 1
    f = lambda x, y, t: 0
    V = lambda x, y: 0
    exact = lambda x, y, t: A*np.cos(kx*x)*np.cos(ky*y)*np.cos(w*t)
    h = [0.1, 0.05, 0.01, 0.005]
    rate = np.zeros(4)
    diff = np.zeros(4)

    # calculate convergence
    for i in range(4):
        dt = 0.1*h[i]
        Nx = int(round(float(Lx)/h[i])) - 2
        Ny = int(round(float(Ly)/h[i])) - 2
        diff[i] = scalar(I, V, f, q, b, Lx, Nx, Ly, Ny, dt, T, exact, user_action='diff', plot=False)
        if i > 0:
            rate[i] = np.log((diff[i]/diff[i-1]))/np.log((float(h[i])/h[i-1]))

    return rate, h

if __name__ == '__main__':
    rate, h = standing_undamped(1, 1., 1.)
    print '  rate   |  dt  '
    for i in range(len(h)):
        print '%.5g   | %g  ' % (rate[i], 0.1*h[i])

