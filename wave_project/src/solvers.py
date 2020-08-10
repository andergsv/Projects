import numpy as np
import sympy as sym
import time as time
from math import exp,sqrt,sin,cos,pi
from scitools.std import *

'''
These functions simulate a 2D wave propagation with neumann boundary conditions
The initial parameters are given by the functions I(x,y) and V(x,y), the 
position and velocity respectively. 
f(x,y,t) is the source term of the equation.
q(x,y) is the function representing the wave speed over the media.
b is the damping factor
Lx and Ly are the limits of the mesh (note: x and y both starf from 0)
Nx and Ny are the number of mesh points for each axis
dt is the time step (it has to be chosen in order to make the result stable)
T is the time of simulation
plot is a flag to set to have a visual plot of the function over time or not
'''



def vectorized(I, V, f, q, b, Lx, Nx, Ly, Ny, dt, T, plot=False):
    # space
    x = np.linspace(0-Lx/Nx, Lx+Lx/Nx, Nx+3)
    dx = x[1] - x[0]
    y = np.linspace(0-Ly/Ny, Ly+Ly/Ny, Ny+3)
    dy = y[1] - y[0]

    #auxiliary variables
    Cx = (dt/dx)**2
    Cy = (dt/dy)**2
    aux1 = 1./(1+0.5*b*dt)
    aux2 = (0.5*b*dt - 1)
    dt2 = dt**2

    # time
    Nt = int(round(float(T)/dt))
    T = Nt*dt
    t = np.linspace(0, T, Nt+1)

    #function matrices
    u = np.zeros([Nx+3, Ny+3])
    u1 = np.zeros([Nx+3, Ny+3])
    u2 = np.zeros([Nx+3, Ny+3])

    #index arrays
    Ix = range(1, u.shape[0]-1)
    Iy = range(1, u.shape[1]-1)
    It = range(0, t.shape[0])
    
    #matrices for parameters q,f,V needed for vectorization
    q_ = np.zeros([Nx+3, Ny+3])
    for i in range(0, u.shape[0]):
        for j in range(0, u.shape[1]):
            q_[i,j]=q(x[i],y[j]);
    V_ = np.zeros([Nx+3, Ny+3])
    V_[:,:]=V(x[:,newaxis],y[newaxis,:])
    f_ = np.zeros([Nx+3, Ny+3])
    for i in range(0, u.shape[0]):
        for j in range(0, u.shape[1]):
            f_[i,j]=f(x[i],y[j],t[0]);
    
    #initial conditions
    u1[1:-1,1:-1]=I(x[1:-1,newaxis],y[newaxis,1:-1])
    ghost_cells(u1,Ix[0],Iy[0],Ix[Nx],Iy[Ny])

    #first time step: 
    u[1:-1, 1:-1] = u1[1:-1, 1:-1]  + \
                          0.5*Cx*(0.5*(q_[1:-1,1:-1] + q_[2:,1:-1])*(u1[2:,1:-1] - u1[1:-1, 1:-1]) - \
                                  0.5*(q_[1:-1,1:-1] + q_[:-2,1:-1])*(u1[1:-1, 1:-1] - u1[:-2, 1:-1]) )+ \
                          0.5*Cy*(0.5*(q_[1:-1,1:-1] + q_[1:-1,2:])*(u1[1:-1,2:] - u1[1:-1, 1:-1]) - \
                                  0.5*(q_[1:-1,1:-1] + q_[1:-1,:-2])*(u1[1:-1, 1:-1] - u1[1:-1, :-2]) )+ \
                          0.5 *dt2*f_[1:-1,1:-1]+ \
                          dt*V_[1:-1,1:-1]

    ghost_cells(u,Ix[0],Iy[0],Ix[Nx],Iy[Ny])
    u2[:,:], u1[:,:] = u1, u

    if plot == True:
        mesh(x[1:-1], y[1:-1], u[1:-1,1:-1], title='t=%g' % t[0], zlim=[-1,1])
        time.sleep(2)

    #main loop 
    for n in It[1:-1]:
        #source term matrix f
        for i in range(0, u.shape[0]):
            for j in range(0, u.shape[1]):
                f_[i,j]=f(x[i],y[j],t[n]);
        #inner points
        u[1:-1, 1:-1] = aux1 * ( \
                          aux2*u2[1:-1, 1:-1]  + 2*u1[1:-1, 1:-1]  + \
                          Cx*(0.5*(q_[1:-1,1:-1] + q_[2:,1:-1])*(u1[2:,1:-1] - u1[1:-1, 1:-1]) - \
                              0.5*(q_[1:-1,1:-1] + q_[:-2,1:-1])*(u1[1:-1, 1:-1] - u1[:-2, 1:-1]) )+ \
                          Cy*(0.5*(q_[1:-1,1:-1] + q_[1:-1,2:])*(u1[1:-1,2:] - u1[1:-1, 1:-1]) - \
                              0.5*(q_[1:-1,1:-1] + q_[1:-1,:-2])*(u1[1:-1, 1:-1] - u1[1:-1, :-2]) )+ \
                          dt2*f_[1:-1,1:-1] )
   
        ghost_cells(u,Ix[0],Iy[0],Ix[Nx],Iy[Ny])
        u2[:,:], u1[:,:] = u1, u

        if plot==True:
            mesh(x[1:-1], y[1:-1], u[1:-1,1:-1], title='t=%g' % t[n], zlim=[-1,1])
            time.sleep(0.01)





def scalar(I, V, f, q, b, Lx, Nx, Ly, Ny, dt, T, ue=None, user_action=False,  plot=False):
    # space
    x = np.linspace(0-Lx/Nx, Lx+Lx/Nx, Nx+3)
    dx = x[1] - x[0]
    y = np.linspace(0-Ly/Ny, Ly+Ly/Ny, Ny+3)
    dy = y[1] - y[0]

    #auxiliary variables
    Cx = (dt/dx)**2
    Cy = (dt/dy)**2
    aux1 = 1./(1+0.5*b*dt)
    aux2 = (0.5*b*dt - 1)
    dt2 = dt**2

    # time
    Nt = int(round(float(T)/dt))
    T = Nt*dt
    t = np.linspace(0, T, Nt+1)

    #function matrices
    u = np.zeros([Nx+3, Ny+3])
    u1 = np.zeros([Nx+3, Ny+3])
    u2 = np.zeros([Nx+3, Ny+3])

    #index arrays
    Ix = range(1, u.shape[0]-1)
    Iy = range(1, u.shape[1]-1)
    It = range(0, t.shape[0])

    #initial conditions
    for i in Ix:
        for j in Iy:
            u1[i][j] = I(x[i], y[j])
    ghost_cells(u1,Ix[0],Iy[0],Ix[Nx],Iy[Ny])


    
    # if actions is chosen
    if user_action == 'diff':
        u_e = np.zeros([Nx+1, Ny+1])
        diff = 0
        for i in Ix:
            for j in Iy:
                u_e[i-1][j-1] = ue(x[i],y[j],t[0])
        uemax = u_e[:].max()
        umax = (u1[1:-1, 1:-1]).max()
        diff = abs(uemax - umax)
    
    if user_action == 'diff2':
        diff2 = 0
        uemax = ue(x[i],y[j],t[0])
        umax = (u1[1:-1, 1:-1]).max()
        diff2 = abs(uemax - umax)
    
    if user_action == 'plug':
        maximum = []
        minimum = []
        maximum.append(u1[:].max())
        minimum.append(u1[:].min())
 
    #first time step:
    for i in Ix:
        for j in Iy:
            qx = 0.5*(q(x[i],y[j]) + q(x[i+1],y[j]))*(u1[i+1][j] - u1[i][j]) - \
                 0.5*(q(x[i],y[j]) + q(x[i-1],y[j]))*(u1[i][j] - u1[i-1][j])
            qy = 0.5*(q(x[i],y[j]) + q(x[i],y[j+1]))*(u1[i][j+1] - u1[i][j]) - \
                 0.5*(q(x[i],y[j]) + q(x[i],y[j-1]))*(u1[i][j] - u1[i][j-1])
            u[i][j] =  u1[i][j] + 0.5*Cx* qx + 0.5*Cy*qy + 0.5*dt**2*f(x[i],y[j], t[0])+V(x[i],y[j])*dt  
    
    ghost_cells(u,Ix[0],Iy[0],Ix[Nx],Iy[Ny])

    if user_action == 'plug':
        maximum.append(u[:].max())
        minimum.append(u[:].min())

    if user_action == 'diff':
        for i in Ix:
            for j in Iy:
                u_e[i-1][j-1] = ue(x[i],y[j],t[0])
        uemax = u_e[:].max()
        umax = (u[1:-1, 1:-1]).max()
        new_diff = abs(uemax - umax)
        if new_diff > diff:
            diff = new_diff 
    
    if user_action == 'diff2':
        uemax = ue(x[i],y[j],t[0])
        umax = (u[1:-1, 1:-1]).max()
        new_diff2 = abs(uemax - umax)
        if new_diff2 > diff2:
            diff2 = new_diff2         

    u2[:], u1[:] = u1, u

    if plot==True:
        mesh(x[1:-1], y[1:-1], u[1:-1,1:-1], title='t=%g' % t[0], zlim=[-1,1])
        time.sleep(2)

    #main loop
    for n in It[1:-1]:
        #inner points
        for i in Ix:
            for j in Iy:
                qx = 0.5*(q(x[i],y[j]) + q(x[i+1],y[j]))*(u1[i+1][j] - u1[i][j]) - \
                     0.5*(q(x[i],y[j]) + q(x[i-1],y[j]))*(u1[i][j] - u1[i-1][j])
                qy = 0.5*(q(x[i],y[j]) + q(x[i],y[j+1]))*(u1[i][j+1] - u1[i][j]) - \
                     0.5*(q(x[i],y[j]) + q(x[i],y[j-1]))*(u1[i][j] - u1[i][j-1])
                u[i][j] = aux1 * ( \
                          aux2*u2[i][j] + 2*u1[i][j] + \
                          Cx*qx + Cy*qy + dt**2*f(x[i],y[j], t[n]) )
              
        ghost_cells(u,Ix[0],Iy[0],Ix[Nx],Iy[Ny])

        if user_action == 'plug':
            maximum.append(u[:].max())
            minimum.append(u[:].min())

        if user_action == 'diff':
            for i in Ix:
                for j in Iy:
                    u_e[i-1][j-1] = ue(x[i],y[j],t[n])
            #print u_e[:]
            #print u[1:-1, 1:-1]
            uemax = u_e[:].max()
            umax = (u[1:-1, 1:-1]).max()
            
            new_diff = abs(uemax - umax)
            if new_diff > diff:
                diff = new_diff
        
        if user_action == 'diff2':
            #print u_e[:]
            #print u[1:-1, 1:-1]
            uemax = ue(x[i],y[j],t[n])
            umax = (u[1:-1, 1:-1]).max()
            
            new_diff2 = abs(uemax - umax)
            if new_diff2 > diff2:
                diff2 = new_diff2
                
        u2[:], u1[:] = u1, u
        
        if plot==True:
            mesh(x[1:-1], y[1:-1], u[1:-1,1:-1], title='t=%g' % t[n], zlim=[-1,1])
            time.sleep(0.01)

    if user_action == 'plug':
        return maximum, minimum
    if user_action == 'diff':
        return diff
    if user_action == 'diff2':
        return diff2


def ghost_cells(v,ix,iy,ixl,iyl):
    #ghost cells for borders  
    v[ix-1,1:-1] = v[ix+1,1:-1] 
    v[ixl+1,1:-1] = v[ixl-1,1:-1] 
    v[1:-1,iy-1] = v[1:-1,iy+1]
    v[1:-1,iyl+1] = v[1:-1,iyl-1]
 
    #ghost cells for corners
    v[ix-1][iy-1] = v[ix+1][iy+1]
    v[ix-1][iyl+1] = v[ix+1][iyl-1]
    v[ixl+1][iy-1] = v[ixl-1][iy+1]
    v[ixl+1][iyl+1] = v[ixl-1][iyl-1]

