from solvers import *
import numpy as np

def plug(direction, c, b, Lx, Ly, T):

    Nx = 20
    Ny = 20
    f = lambda x, y, t: 0
    q = lambda x, y: c**2
    V = lambda x, y: 0
    amplitude = 2.0

    if direction == 'x':
        I = lambda x, y: amplitude if x < 1 else 0
    else:
        I = lambda x, y: amplitude if y < 1 else 0
    
    # control value of dt in order
    # to make sure c*dt/dx=1
    dt = float(Lx)/(c*(Nx+2))

    maximum, minimum = scalar(I, V, f, q, b, Lx, Nx, Ly, Ny, dt, T, user_action='plug' , plot=True)
    
    if max(maximum) == amplitude and min(maximum) == amplitude/2:
        print 'test ok for maximum values!'
    else:
        print 'test failed for minimum values'
    if any(minimum) != 0.0:
        print 'test failed for minimum values'
    else:
        print 'test ok for minimum values!'
            
    
if __name__ == '__main__':
    plug(direction='x',c=1, b=0, Lx=10, Ly=10, T=10)
    plug(direction='y',c=1, b=0, Lx=10, Ly=10, T=10)

