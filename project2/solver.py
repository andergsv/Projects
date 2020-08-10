from dolfin import *
import numpy

global error
def solver(N, dim, Pdeg, dt, u0, f, rho, alpha, exercise):
    set_log_active(False)
    '''
    Neumann bc already standard in fenics
    '''
    nx = ny = nz = N
    # Create mesh and define function space
    if dim==1:
        mesh = UnitIntervalMesh(nx)
    elif dim ==2:
        mesh = UnitSquareMesh(nx, ny)
    elif dim ==3:
        mesh = UnitCubeMesh(nx, ny, nz) 
    
    # Defineing a discrete function space
    V = FunctionSpace(mesh, 'CG', Pdeg) #Continuous Galerkin
    
    # Initial condition
    u_1 = interpolate(u0, V)

    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    
    a = u*v*dx + dt/rho*inner(alpha(u_1)*nabla_grad(u), nabla_grad(v))*dx
    L = (u_1 + dt/rho*f)*v*dx

    # Compute solution
    u = Function(V)   # the unknown at a new time level
    T = 1             # total simulation time
 
    t = dt    
    while t <= T:
        u0.t = t
        f.t = t
        solve(a == L, u)

        # Verify
        u_e = interpolate(u0, V)
        if exercise == 'd':
            maxdiff = numpy.abs(u_e.vector().array() - u.vector().array()).max()
            print('Max error, t=%.2f: %-10.3f' % (t, maxdiff))
        if exercise == 'e':
            if abs(t-0.5) < 0.001:
                e = u_e.vector().array() - u.vector().array()
                E = numpy.sqrt(numpy.sum(e**2)/u.vector().array().size)
                print' h=%.5f, E/h=%.5f, t=%0.2f' %(numpy.sqrt(dt),E, t)  
                
        if exercise == 'f':
            import matplotlib.pyplot as plt
            x = numpy.linspace(0,1,u.vector().array().size)
            plt.plot(x, u_e.vector().array()[::-1],'--', label='exact')
            plt.plot(x, u.vector().array()[::-1],'o',label='numerical')
            plt.xlabel('x')
            plt.ylabel('u(x)')
            plt.legend(loc='upper left')
            plt.title('t=%s' %t)
            plt.savefig('tex/ft_'+str(int(t*100))+'.png')
            #plt.axis([x[0], x[-1], 0 ,0.2])
            plt.show()
        
        if exercise == 'h':
            e = u_e.vector().array() - u.vector().array()
            E = numpy.sqrt(numpy.sum(e**2)/u.vector().array().size).max()
            return E
        if exercise == 'i':
            from time import sleep
            umax = u.vector().array().max()
            umin = u.vector().array().min()
            if umax - umin < 0.00001:
                wiz = plot(u)#, interactive = True)
                wiz.write_png('t_'+str(t))
                break
                
            elif t == dt:
                wiz = plot(u)#, interactive = True)
                wiz.write_png('t_'+str(t))
            
            else:
                plot(u).set_min_max(0.3,1)
        t += dt
        u_1.assign(u)
        
        
def d():
    u0 = Expression('3.0',t=0)
    alpha = lambda u: 5 #Expression('5', u)
    f = Constant('0')
    rho = 1000
    dim = [1,2,3]
    Pdeg = [1,2]
    N = 5
    dt = 0.5
    for i in Pdeg:
        print '-------P%s element-------' %i
        for j in dim:
            print '----------%s D-----------' %j
            solver(N, j, i, dt, u0, f, rho, alpha, 'd')
            
def e():
    u0 = Expression('exp(-pi*pi*t)*cos(pi*x[0])',t=0)
    alpha = lambda u: 1
    f = Constant('0')
    rho = 1
    dim = 2
    Pdeg = 1 
    dt = numpy.zeros(6)
    dt[0] = 0.5
    for i in range(len(dt)-1):
        dt[i+1] = dt[i]/2.
    N = [int(round(1./numpy.sqrt(j))) for j in dt] 
    for j in range(len(N)):
        solver(N[j], dim, Pdeg, dt[j], u0, f, rho, alpha, 'e')

def f():
    u0 = Expression('t*pow(x[0],2)*(0.5 - x[0]/3.)',t=0)
    alpha = lambda u: 1 + u**2
    rho = 1
    f = Expression('-rho*pow(x[0],3)/3 + rho*pow(x[0],2)/2 + 8*pow(t,3)*pow(x[0],7)/9 - \
                    28*pow(t,3)*pow(x[0],6)/9 + 7*pow(t,3)*pow(x[0],5)/2 - \
                    5*pow(t,3)*pow(x[0],4)/4 + 2*t*x[0] - t', rho=rho, t=0)
    N = 10
    dim = 1
    Pdeg = 1 
    dt = 0.2
    solver(N, dim, Pdeg, dt, u0, f, rho, alpha, 'f')

def h():
    u0 = Expression('t*pow(x[0],2)*(0.5 - x[0]/3.)',t=0)
    alpha = lambda u: 1 + u**2
    rho = 1
    dim = 1
    Pdeg = 1 
    dt = [0.5, 0.04, 0.0024, 0.0003, 0.00002, 0.000002]
    N = [int(round(1./numpy.sqrt(k))) for k in dt] 
    er = []
    for j in range(len(N)):
        error = 0
        f = Expression('rho*pow(x[0],2)*(-2*x[0] + 3)/6. - \
                        (-12*t*x[0] + 3*t*(-2*x[0] + 3))* \
                        (pow(x[0],4)*pow((- dt + t),2)*pow((-2*x[0] + 3),2) + \
                        36)/324.- (-6*t*pow(x[0],2) + 6*t*x[0]*(-2*x[0] + 3))* \
                        (36*pow(x[0],4)*pow((- dt + t),2)*(2*x[0] - 3) + \
                        36*pow(x[0],3)*pow((- dt + t),2)* pow((-2*x[0] + 3),2))/5832.', 
                        rho=rho, t=0, dt=dt[j])
        E = solver(N[j], dim, Pdeg, dt[j], u0, f, rho, alpha, 'h')
        if E > error:
            error = E
            
        er.append(error)
        
    for i in range(len(er)-1):
        convrate = numpy.log(abs(er[i+1] - er[i]))/numpy.log(abs(dt[i+1] - dt[i]))
        print 'h=%.6f,   convergence rate:%s' %(dt[i+1], convrate)

def i():
    sigma = 0.5
    beta = 9
    f = Constant('0')
    u0 = Expression('exp(-1./(2*pow(sigma,2))*(pow(x[0],2)+pow(x[1],2)))',t=0, sigma=sigma)
    alpha = lambda u: 1 + beta*u**2
    rho = 1
    dim = 2
    Pdeg = 1 
    dt = 0.005
    N = 10
    solver(N, dim, Pdeg, dt, u0, f, rho, alpha, 'i')    

if __name__ == '__main__':
    #d()
    #e()
    #f()
    h()
    #i()
