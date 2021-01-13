#u_t + f(u)_x = 0, f(u) = u^2/(u^2+a(1-u)^2), i.e., the Buckley-Leverett model

import numpy as np
import matplotlib.pyplot as plt
import argparse

a_buck = 1.0 #The constant a in definition of flux
u_buck = 0.5 #Depends on a, we have f convex in (0,u_buck),concave in (u_buck,1)

def flux(u): #Buckley-Leverett flux
  a = a_buck
  num = u**2;den=u**2+a*(1.-u)**2
  return num/den

def num_flux(ul,ur):
  if ul <= ur:
    z = np.linspace(ul,ur,100)
    return flux(z).min()
  else:
    z = np.linspace(ur,ul,100)
    return flux(z).max()

#h*du/dt = res(u); res(u) = g_{}
def compute_residual(u):
  n = len(u)
  res = np.zeros(n)
  for j in range(0,n-1): #each j computes flux g_{j+1/2}
    flux      = num_flux(u[j],u[j+1])
    res[j]   -= flux
    res[j+1] +=flux
  #last face
  flux      = num_flux(u[n-2],u[n-1])#f_{n-1/2}
  res[n-1] -= flux
  res[0]   += flux
  return res
# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of cells', default=100)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.4)
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)
args = parser.parse_args()
N= args.N
dx = 1./N;cfl = args.cfl;dt = cfl*dx
xmin,xmax = -1,1.

x = np.linspace(xmin,xmax,N)
'''
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
  
def exact_soln(x,t):
  f = np.empty_like(x)
  u_star = 1.2071067811865481
  for i,xx in enumerate(x):
    if xx <= 0.*t:#x<f'(1)t
      f[i] = 1.
    elif xx > 0.*t and xx < u_star*t:
      f[i] = 1.*(xx-u_star*t)/(0.*t-u_star*t)#Not correct, temporary
    elif xx > u_star*t:
      f[i] = 0.
  return f

#u = np.ones(N)
u = initial_condition(x)
t, it = 0.0, 0
Tf = args.Tf

uold = np.empty_like(u)
fig = plt.figure()
ax = fig.add_subplot(111)
line1, = ax.plot(x, u, 'ro')
line2, = ax.plot(x, u, 'b')
ax.set_xlabel('x'); ax.set_ylabel('u')
plt.legend(('Numerical','Exact'))
plt.title('N='+str(N)+', CFL='+str(cfl)+', Scheme= Godunov')
plt.axis([xmin, xmax, u.min()-0.1, u.max()+0.1])
plt.grid(True); plt.draw(); plt.pause(0.5/N)
wait = input("Press enter to continue ")
lam = dt/dx
while t < Tf:
    uold[:] = u
    res = compute_residual(u)
    u = uold + lam*res
    t += dt; it += 1
    line1.set_ydata(u)
    line2.set_ydata(exact_soln(x,t))
    plt.draw(); plt.pause(0.05)
plt.show()

def flux_buck(u):
  a = a_buck
  num = u**2;denom=u**2+a*(1.-u)**2
  return num/denom

