from solver import *
from dolfin import *
from mshr import *
T = 2000
alpha = 0.1 #Amplitude factor: >0 Weakly non-linear, =0 linear (0.3/0.1 in report)

############# Eta Parameters##############
c = lambda a: np.sqrt(6*(1+a)**2/(a**2*(3+2*a))*((1+a)*np.log(1+a)-a))
beta = lambda a: np.sqrt(3*a/(4*(1+0.68*a)))

    
#rec = Rectangle(Point(0,0), Point(20, 10))
#cyl = Circle(Point(15+3/2., 3.5+3/2.), 3/2., 20)
#mesh = generate_mesh(rec - cyl, 23) 
#alpha, eps, dt, mesh, H, kappa, T

width, length = 5, 50
P2 = (length,width)
P3 = (P2[0]+P2[1],40)
L = lbend(50,20,20,50,P2,P3)
mesh = L.mesh()
plot(mesh,interactive=True)
pl = FEMbouss(1.2,1.2, 0.1, mesh, Constant(1), 1e7, T)

#Expression('1 - 0.5*x[0]')

A = 1. 
i = 1
j = 0
kxc = i*np.pi
kyc = j*np.pi
mu = 1
delta = 0
w2 = (kxc**2+kyc**2)/(1+mu**2*(kxc**2+kyc**2)/3.)
B = -float(A)*np.sqrt(w2)/(kxc**2+kyc**2)


x0 = 0
eta00 = Expression('a*4*pow(cosh(B*(x[0]-x0)),2)/pow(cosh(2*B*(x[0]-x0)) + 1,2)/(1+a*pow(tanh(B*(x[0]-x0)),2))', a = alpha, B = beta(alpha), x0 = x0, element=FiniteElement("Lagrange", triangle, 1))
    
      
    
initeta = Expression('A*cos(w*t+delta)*cos(kx*x[0])*cos(ky*x[1])', A = A, delta=delta, t=0, w=np.sqrt(w2), kx=kxc, ky=kyc, degree=1)
initu = Expression('-B*kx*sin(w*t+delta)*sin(kx*x[0])*cos(ky*x[1])', B = B, delta=delta, t=0, w=np.sqrt(w2), kx=kxc, ky=kyc, degree=1)
initv = Expression('-B*ky*sin(w*t+delta)*cos(kx*x[0])*sin(ky*x[1])', B = B, delta=delta, t=0, w=np.sqrt(w2), kx=kxc, ky=kyc,degree=1)
phiinit = Expression('B*sin(w*t+delta)*cos(kx*x[0])*cos(ky*x[1])', B = B, delta=delta, t=0, w=np.sqrt(w2), kx=kxc, ky=kyc, degree=1)

pl.usolver(eta00, initu, initv)
#eta2 = pl.potsolver(initeta, phiinit)

#pl.globouss()
#print errornorm(eta1, eta2, 'l2')
