#u_t + f(u)_x = 0, f(u) = u^2/(u^2+a(1-u)^2), i.e., the Buckley-Leverett model

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import argparse

a_buck = 1.0 #The constant a in definition of flux
u_buck = 0.5 #Depends on a, we have f convex in (0,u_buck),concave in (u_buck,1)
u_s,u_ss = 1./np.sqrt(2.),1.-1./np.sqrt(2.)
def flux(u): #Buckley-Leverett flux
  a = a_buck
  num = u**2;den=u**2+a*(1.-u)**2
  return num/den

def fprime(u):
    a = a_buck
    L = u**2 + a*(1.0-u)**2
    return 2.0*a*u*(1.0-u) / L**2

#Inverse of restriction of f' to [1/sqrt(2),1]
def inv_f_s(v):
    def f(x):
        p = 4.*v*x**4-8.*v*x**3+(8.*v+2.)*x**2+(-4.*v-2.)*x+v
        return p
    return optimize.brentq(f, u_s, 1.,xtol=2e-12, rtol=8.881784197001252e-16)
#Inverse of restriction of f' to [0.,1-1/sqrt(2)]
def inv_f_ss(v):
    def f(x):
        p = 4.*v*x**4-8.*v*x**3+(8.*v+2.)*x**2+(-4.*v-2.)*x+v
        return p
    return optimize.brentq(f, 0., u_ss,xtol=2e-12, rtol=8.881784197001252e-16)

def num_flux(ul,ur):
  return flux(ul)#Since flux is increasing
  '''
  if ul <= ur:
    #z = np.linspace(ul,ur,100)
    return flux(ul)
  else:
    #z = np.linspace(ur,ul,100)
    return flux(ul)
  '''
    




# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of cells', default=100)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.4)
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)
parser.add_argument('-ic', type=str, help = 'Initial condition',default='hatbuck')
parser.add_argument('-bc',type = str,default = 'periodic')
args = parser.parse_args()



N= args.N
#WHAT A bug. I had put dx = 1./N, and that caused the solution to run at twice speed!
xmin,xmax = -1.,1.
dx = (xmax-xmin)/N;cfl = args.cfl;dt = cfl*dx

x = np.linspace(xmin,xmax,N)
def initial_condition(x):
    f = np.empty_like(x)
    length = xmax-xmin
    for i,xx in enumerate(x):
        if xx < xmin + 0.25*length or xx > xmin+0.5*length:
            f[i] = 0.0
        else:
            f[i] = 1.0
    return f
'''
def initial_condition(x):
    f = np.empty_like(x)
    length = xmax-xmin
    for i,xx in enumerate(x):
        if xx < xmin + 0.5*length:
            f[i] = 1.0
        else:
            f[i] = 0.0
    return f
'''
if args.ic=='step':
  #Not exact soln, but has shock and rarefaction at same place
  def exact_soln(x,t):
    f = np.empty_like(x)
    u_star = 1./np.sqrt(2.)#satisfied f(u_star)-f(0)/(u_star-0)=f'(u_star)
    f_u_star = fprime(u_star) #f'(u_star)
    for i,xx in enumerate(x):
      if xx <= 0.*t:#x<f'(1)t
        f[i] = 1.
      elif xx > 0.*t and xx < f_u_star*t:
        f[i] = inv_f_s(xx/t)
      elif xx > f_u_star*t:
        f[i] = 0.
    return f
else:
  #Not exact soln, but has shock and rarefaction at same place
  def exact_soln(x,t):
    f = np.empty_like(x)
    u_s,u_ss = 1./np.sqrt(2.),1.-1./np.sqrt(2.)
    #f((u_s)-f(0))/(u_s-0)=f'(u_s),(f(u_ss)-f(1))/(u_ss-1)=f'(u_ss)
    f_u_s,f_u_ss = fprime(u_s),fprime(u_ss) #f'(u_star)
    for i,xx in enumerate(x):
      if xx <= -0.5:
        f[i] = 0.
      elif xx>-0.5 and xx <= -0.5+f_u_ss*t:
        f[i] = inv_f_ss((xx-(-0.5))/t)
      elif xx>=-0.5+f_u_ss*t and xx <= 0.:
        f[i] = 1.
      elif xx>=0 and xx>= -0.5+f_u_ss*t and xx<=f_u_s*t:
        f[i] = inv_f_s(xx/t)
      else:
        f[i] = 0.
    return f

#h*du/dt = res(u); res(u) = g_{}
def compute_residual(u):
  n = len(u)
  res = np.zeros(n)
  flux      = num_flux(u[0],u[1])#f_{1/2}
  res[0]   -= flux
  res[1]   += flux
  for j in range(1,n-1): #each j computes flux g_{j+1/2}
    flux      = num_flux(u[j],u[j+1])
    res[j]   -= flux
    res[j+1] +=flux
  #last face
  flux      = num_flux(u[n-2],u[n-1])#f_{n-1/2}
  res[n-1] -= flux
  res[0]   += flux
  return res

#u = np.ones(N)
u = exact_soln(x,0.)
t, it = 0.0, 0
Tf = args.Tf

uold = np.empty_like(u)
fig = plt.figure()
ax = fig.add_subplot(111)
line1, = ax.plot(x, u, 'ro')
if args.bc == 'dirichlet':
  line2, = ax.plot(x, u, 'b')
ax.set_xlabel('x'); ax.set_ylabel('u')
plt.legend(('Numerical','Exact'))
plt.axis([xmin, xmax, u.min()-0.1, u.max()+0.1])
plt.grid(True); plt.draw(); plt.pause(0.1/N)
wait = input("Press enter to continue ")
lam = dt/dx
while t < Tf:
    uold[:] = u
    exact = exact_soln(x,t+dt)
    if args.bc == 'dirichlet':
      u[0] = exact[0]
      u[-1] = u[-2]
    res = compute_residual(u)
    u = uold + lam*res
    t += dt; it += 1
    plt.title('N='+str(N)+', CFL='+str(cfl)+', Scheme= Godunov, '+'t='+str(round(t,3)))
    line1.set_ydata(u)
    if args.bc == 'dirichlet':
      line2.set_ydata(exact)
    plt.draw(); plt.pause(0.1/N)
plt.show()
