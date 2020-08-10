from dolfin import *
import numpy as np
import copy as cp
set_log_active(False)
class lbend:
    def __init__(self, Nx1, Ny1, Nx2, Ny2, Point2, Point3):
        self.Nx1, self.Ny1, self.Nx2, self.Ny2, self.Point1, self.Point2, self.Point3 = Nx1, Ny1, Nx2, Ny2, (0,0), Point2, Point3

    def mesh(self):
        '''
        Nx1, Ny1: Elements first rectangle
        Nx2, Ny2: Elements third rectangle

        Second rectangle uses hy1 and hx2

        Point1 to 2 is the diagonal in first rectangle
        Point2 to 3 is the diagonal in third rectangle

                                       __________x Point 3
                                      |          |
                                      |   Rec. 3 |
                               Point 2|          | 
                ______________________x__________|
                |                     |          |
                |      Rec. 1         |  Rec. 2  |
                |                     |          |
        Point 1 x_____________________|__________|
           (0,0)
        '''

        self.Nx1 = self.Nx1 + 1
        self.Ny1 = self.Ny1 + 1
        self.Nx2 = self.Nx2 + 1
        self.Ny2 = self.Ny2 + 1

        if len(self.Point2) != 2 or len(self.Point1) != 2 or len(self.Point3) != 2:
            print '(x,y) Point needs to be a list with two arguments, x and y'

        xnodes1 = np.linspace(self.Point1[0], self.Point2[0], self.Nx1)
        ynodes1 = np.linspace(self.Point1[1], self.Point2[1], self.Ny1)

        xnodes3 = np.linspace(self.Point2[0], self.Point3[0], self.Nx2)
        ynodes3 = np.linspace(self.Point2[1], self.Point3[1], self.Ny2)

        xnodes2 = xnodes3 # np.linspace(Point2[0], Point3[0], Nx2)
        ynodes2 = ynodes1 # np.linspace(Point1[1], Point2[2], Ny1)


        mergedxlist = list(xnodes1) + list(xnodes2[1:])

        output = open('lmesh.xml', 'w')
        output.write('''<?xml version="1.0" encoding="UTF-8"?>

<dolfin xmlns:dolfin="http://www.fenicsproject.org">
  <mesh celltype="triangle" dim="2">\n''')
        output.write('    <vertices size="%s">\n' %int(len(xnodes1)*len(ynodes1)+len(xnodes2)*len(ynodes2)+len(xnodes3)*len(ynodes3)-len(ynodes1)-len(xnodes3)))
        vertex_index = 0
        for y in ynodes1:
            for x in mergedxlist:
                output.write('      <vertex index="%s" x="%.16e" y="%.16e" z="0.0000000000000000e+00"/>\n' %(int(vertex_index),x,y))
                vertex_index +=1

        for y in ynodes3[1:]:
            for x in xnodes3:
                output.write('      <vertex index="%s" x="%.16e" y="%.16e" z="0.0000000000000000e+00"/>\n' %(int(vertex_index),x,y))
                vertex_index +=1

        numelements = ((self.Nx1-1)*(self.Ny1-1) + (self.Nx2-1)*(self.Ny2-1) + (self.Ny1-1)*(self.Nx2-1))*2
        output.write('\n    </vertices>\n')	
        output.write('    <cells size="%s">\n' %int(numelements))

        v0 = 0
        a = 0
        v1 = len(mergedxlist)+1
        indx = 0

        for j in range(len(ynodes1)-1):
            for i in range(len(mergedxlist)-1):
                output.write('      <triangle index="%s" v0="%s" v1="%s" v2="%s"/>\n' %(indx,v0,v1,v0+1))
                indx += 1
                output.write('      <triangle index="%s" v0="%s" v1="%s" v2="%s"/>\n' %(indx,v0,v1,v1-1))
                indx += 1

                v0+=1
                v1+=1

            v0 += 1
            v1 += 1

        rec3 = v1-self.Nx2-1
        for k in range(len(ynodes3)-1):
            for r in range(len(xnodes3)-1):
                output.write('      <triangle index="%s" v0="%s" v1="%s" v2="%s"/>\n' %(indx,rec3,v1,rec3+1))
                indx += 1
                output.write('      <triangle index="%s" v0="%s" v1="%s" v2="%s"/>\n' %(indx,rec3,v1,v1-1))
                indx += 1

                rec3+=1
                v1+=1
			
            rec3 += 1
            v1 += 1
        output.write('''    </cells>
  </mesh>
</dolfin>
    ''')
        output.close()
        return Mesh('lmesh.xml')


    def refine(self):
        self.Nx1, self.Ny1, self.Nx2, self.Ny2 = 2*(self.Nx1-1), 2*(self.Ny1-1), 2*(self.Nx2-1), 2*(self.Ny2-1)
        return self.mesh()


class FEMbouss:
    def __init__(self, alpha, eps, dt, mesh, H, kappa, T):
        self.alpha = alpha                                  #amplitudefaktor
        self.eps = eps                                      #Langbolgeutvikler
        self.dt = dt
        self.mesh = mesh
        self.H = H                                       #Depth function/(next step: matrix)
        self.kappa = kappa                               #Nitsche constant
        self.T = T                                       #End time
    
    def potsolver(self, initial_eta, initial_phi):   #Constant depth, irrotational
        print 'returned values (eta, u):'
        retval = raw_input()
        print retval
        V = FunctionSpace(self.mesh, 'CG', 1)
        
        eta0 = project(initial_eta, V)
        phi0 = project(initial_phi, V)

        phi = TrialFunction(V)
        eta = TrialFunction(V)
        Ni = TestFunction(V)

        phi_ = Function(V)
        eta_ = Function(V)

        F1 = Constant(1/self.dt)*(phi-phi0) * Ni * dx +  eta0 * Ni * dx \
            + Constant(self.alpha*0.5)*dot(grad(phi0),grad(phi))* Ni*dx \
            + Constant(self.eps/(self.dt*3)) * dot(grad(Ni),grad((phi-phi0)))*dx 
        F2 = Constant(1/self.dt)*(eta-eta0) * Ni * dx \
            - inner(grad(phi0), grad(Ni)) * dx \
            - Constant(self.alpha)*eta*dot(grad(phi0),grad(Ni))*dx
            
        a1 = assemble(lhs(F1))
        a2 = assemble(lhs(F2))
        
        t = self.dt
        while t< self.T:
            b1 = assemble(rhs(F1))
            solve(a1, phi_.vector(), b1)
            phi0.assign(phi_)
            
            b2 = assemble(rhs(F2))
            solve(a2, eta_.vector(), b2)
            eta0.assign(eta_)
            plot(eta_, rescale=False)
                
            t +=self.dt
        return eta0
        
    def usolver(self, initial_eta, initial_u, initial_v):
        self.initial_eta = initial_eta
        self.initial_u = initial_u
        self.initial_v = initial_v
        
        '''
        eta as C++ Expression. Further work: implement as matrix
        u and v as expressions
        '''    
        V = VectorFunctionSpace(self.mesh, 'CG', 2)
        Q = FunctionSpace(self.mesh, 'CG', 1) 
        
        v = TestFunction(V)
        q = TestFunction(Q) 

        u = TrialFunction(V)
        eta = TrialFunction(Q)

        n = FacetNormal(self.mesh)
        
        eta0 = project(self.initial_eta,Q)
        u0 = project(Expression(('u','v'), u=self.initial_u, v=self.initial_v, \
                                 degree=2), V)
        h = project(self.H,Q)
        
        self.eta_ = Function(Q)
        self.u_ = Function(V)
        
        F1 = Constant(1/self.dt)*dot((u-u0),v)*dx + dot(grad(eta0),v)*dx \
            + Constant(self.alpha)*inner((dot(u0,grad(u)) \
            + dot(u,grad(u0)))/Constant(2),v)*dx \
            - Constant(self.eps/(self.dt*6))*inner(div((u-u0)),div(h*h*v))*dx \
            + Constant(0.5*self.eps/self.dt)*inner(div(h*(u-u0)),div(h*v))*dx
        F1 +=  Constant(self.eps/(self.dt*6))*inner(div((u-u0)),h*h*dot(v,n))*ds \
            + Constant(self.eps/(self.dt*6))*inner(div(h*h*v),dot((u-u0),n))*ds
            - Constant(0.5*self.eps/self.dt)*inner(div(h*(u-u0)),div(h*v))*dx
            -
            +Constant(self.kappa*self.eps)*h*h/Constant(3*self.dt)*dot(v,n) \
            *(dot(u,n)-dot(u0,n))*ds \
            
        F2 = Constant(1/self.dt)*(eta-eta0)*q*dx \
            + Constant(self.alpha)*(dot(grad(eta),u0) + eta*div(u0))*q*dx \
            - dot(u0,grad(h*q))*dx + dot(grad(h),u0)*q*dx
        
        a1 = assemble(lhs(F1))
        a2 = assemble(lhs(F2))

        self.t = self.dt
        while self.t<self.T:
            b1 = assemble(rhs(F1))
            solve(a1, self.u_.vector(), b1)
            u0.assign(self.u_)
            #plot(u_)
            b2 = assemble(rhs(F2))
            solve(a2, self.eta_.vector(), b2)  
            eta0.assign(self.eta_)
            
            plot(self.eta_, rescale=False) 
                      
            self.t += self.dt

              
      
            

        

