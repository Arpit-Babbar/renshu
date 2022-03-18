#u_t + f(u)_x = 0, f(u) = u^2/(u^2+a(1-u)^2), i.e., the Buckley-Leverett model

#Need to change it to compute exact solution upto Tf at the beginning for easy plotting.

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import integrate
from numpy import sqrt
import argparse

a_buck = 1.0 # The constant a in definition of flux
u_buck = 0.5 # Depends on a, we have f convex in (0,u_buck),concave in (u_buck,1)

# For exact solution, we need u_s,u_ss from convex hull reconstruction s.t.
# (f(u_s)-f(0))/(u_s-0)=f'(u_s), (f(u_ss)-f(1))/(u_ss-1)=f'(u_ss)
u_s,u_ss = 1.0/np.sqrt(2.0),1.0-1.0/np.sqrt(2.0)
def flux(u): #Buckley-Leverett flux
  a = a_buck
  num = u**2;den=u**2+a*(1.-u)**2
  return num/den

def fprime(u): #derivative of flux
    a = a_buck
    L = u**2 + a*(1.0-u)**2
    return 2.0*a*u*(1.0-u) / L**2

# Inverse of f' restricted to [0.5,1], the interval that contains [u_s,1]
# This is used to compute rarefaction solutions
def inv_f_s(v):
    # Inverse of f' at v equals root of this polynomial in [0.5,1]
    def p(x):
        value = 4.*v*x**4-8.*v*x**3+(8.*v+2.)*x**2+(-4.*v-2.)*x+v
        return value
    if 0<=v and v<= 2.:
        output = optimize.brentq(p, 0.5, 1.)# Gives root of polynomial
    else:
        print('Inverse of f\' cannot be computed for v = ',v)
        print('Please give values in [0,2] only')
    return output

# Inverse of f' restricted to [0,0.5], the interval that contains [0,u_ss]
def inv_f_ss(v):
    # Inverse of f' at v equals root of this polynomial in [0,0.5]
    def p(x):
        value = 4.*v*x**4-8.*v*x**3+(8.*v+2.)*x**2+(-4.*v-2.)*x+v
        return value
    if 0<=v and v<= 2.:
        output = optimize.brentq(p, 0., 0.5)# Gives root of polynomial
    else:
        print('Inverse of f\' cannot be computed for v = ',v)
        print('Please give values in [0,2] only')
    return output

# Godunov flux, defined as max/min b/w u_l,u_r computed using monotonicity
def num_flux(ul,ur):
  return flux(ul)

# To compute the exact solution, we need all shock locations for all time
# The shock location between the two rarefaction waves needs to be computed
# by solving the ode s'(t)=(f_l-f_r)/(u_l-u_r) for the particular time
def update_shock(shock,t,dt):
    def rh(shock,t):
      u_l = inv_f_ss((shock+0.5)/t)
      u_r = inv_f_s(shock/t)
      f_l = flux(u_l)
      f_r = flux(u_r)
      dsdt = (f_l-f_r)/(u_l-u_r)
      return dsdt
    # This is the time where rarefaction characteristics intersect
    if t>=1./(2.*fprime(u_ss)):
      time = [t,t+dt]
      output = integrate.odeint(rh,shock,time,rtol = 1e-5)
      return output[1]
    else:
      return 0.0

def exact_soln_step(x,t):
  f = np.empty_like(x)
  u_star = 1./np.sqrt(2.)# satisfied f(u_star)-f(0)/(u_star-0)=f'(u_star)
  f_u_star = fprime(u_star) # f'(u_star)
  for i,xx in enumerate(x):
    if xx <= 0.*t:# x<f'(1)t
      f[i] = 1.
    elif xx > 0.*t and xx < f_u_star*t:
      f[i] = inv_f_s(xx/t)
    elif xx > f_u_star*t:
      f[i] = 0.
  return f

def exact_soln_hatbuck(x,t,shock):
  f = np.empty_like(x)
  f_u_s,f_u_ss = fprime(u_s),fprime(u_ss) # f'(u_s),f'(u_ss)
  for i,xx in enumerate(x):
    if xx <= -0.5:
      f[i] = 0.0
    elif xx>-0.5 and xx <= -0.5+f_u_ss*t:
      f[i] = inv_f_ss((xx-(-0.5))/t)
    elif xx>=-0.5+f_u_ss*t and xx <= 0.0:
      f[i] = 1.0
    elif xx>=0.0 and xx <=f_u_s*t:
      if xx>= shock:
        f[i] = inv_f_s(xx/t)
      else:
        f[i] = inv_f_ss((xx-(-0.5))/t)
    else:
      f[i] = 0.0
  return f

# h*du/dt = res(u); res(u) = -(g_{j+1/2}-g_{j-1/2})
def compute_residual(u):
  n         = len(u)
  res       = np.zeros(n)
  flux      = num_flux(u[0],u[1])# f_{1/2}
  res[0]   -= flux
  res[1]   += flux
  for j in range(1,n-1): # each j computes flux g_{j+1/2}
    flux      = num_flux(u[j],u[j+1])
    res[j]   -= flux
    res[j+1] += flux
  # last face
  flux      = num_flux(u[n-1],u[0])#f_{n-1/2}
  res[n-1] -= flux
  res[0]   += flux
  return res

def solve(xmin,xmax,N,cfl,Tf,bc):
  dx,dt = (xmax-xmin)/N, cfl/N
  lam = dt/dx

  x = np.linspace(xmin,xmax,N)

  u = exact_soln(x,0.,0.)
  t, it = 0.0, 0

  uold = np.empty_like(u)
  fig = plt.figure()
  ax = fig.add_subplot(111)
  line1, = ax.plot(x, u, 'ro', markersize=3)
  if bc == 'dirichlet' or bc == 'periodic_exact':
    line2, = ax.plot(x, u, 'b')
  ax.set_xlabel('x'); ax.set_ylabel('u')
  plt.legend(('Numerical','Exact'))
  plt.axis([xmin, xmax, u.min()-0.1, u.max()+0.1])
  plt.grid(True); plt.draw(); plt.pause(0.1/N)
  wait = input("Press enter to continue ")
  shock = 0.
  # Initial value of shock
  while t < Tf:
      uold[:] = u
      if (t+dt)>Tf:
        dt = Tf-t
      shock = update_shock(shock,t,dt)# Gives shock for solution at time t+dt
      exact = exact_soln(x,t+dt,shock)
      if bc == 'dirichlet':
        u[0] = exact[0]
      res = compute_residual(u)
      u = uold + lam*res
      t += dt; it += 1
      plt.title('N='+str(N)+', CFL='+str(cfl)+', Scheme= Godunov, '+'t='+str(round(t,3)))
      if it%args.plot_freq == 0 or np.abs(Tf-t) < 1.0e-13:
        line1.set_ydata(u)
        if bc == 'dirichlet' or bc == 'periodic_exact':
          line2.set_ydata(exact)
        plt.draw(); plt.pause(0.1/N)
  print('The Linfty error is ',np.max(exact-u))
  print('The L1 error is ',      dx*np.linalg.norm((exact-u),1))
  print('The L2 error is ',sqrt(dx)*np.linalg.norm((exact-u),2))
  plt.show()

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of cells', default=100)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.4)
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)
parser.add_argument('-ic', type=str, help = 'Initial condition. Choices are  \
                                             hatbuck,step',default='hatbuck')
parser.add_argument('-bc',type = str,help = 'Choices are periodic,dirichlet. \
                    periodic_exact does periodic and plots Dirichlet exact soln', \
                                            default = 'periodic')
parser.add_argument('-plot_freq', type=int, help = 'Plot Frequency')
args = parser.parse_args()

if args.ic=='step':
  def exact_soln(x,t,shock):
    return exact_soln_step(x,t)
elif args.ic=='hatbuck':
  def exact_soln(x,t,shock):
    return exact_soln_hatbuck(x,t,shock)
else:
  print('Incorrect initial condition')
  assert(False)
xmin,xmax=-1.0,1.0
solve(xmin,xmax,args.N,args.cfl,args.Tf,args.bc)
