"""
Solve u_t + u_x = 0 with periodic bc
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse




def update_lw(nu, u):
    unew = np.empty_like(u)
    unew[0] = u[0] - 0.5*nu*(u[1]-u[-2]) + 0.5*nu**2*(u[-2]-2*u[0]+u[1])
    unew[1:-1] = u[1:-1] - 0.5*nu*(u[2:] - u[0:-2]) \
                 + 0.5*nu**2*(u[0:-2] - 2*u[1:-1] + u[2:])
    unew[-1] = unew[0]
    return unew
def update_back_lw(nu,u):
  unew = np.empty_like(u)
  unew[0] = u[0]  - nu*(1.5*u[0]-2.0*u[-2]+0.5*u[-3]) \
            +0.5*nu**2*(u[-3] - 2.0*u[-2] + u[0]) 
  unew[1] =u[1]  - nu*(1.5*u[1]-2.0*u[0]+0.5*u[-2]) \
            +0.5*nu**2*(u[-2] - 2.0*u[0] + u[1])
  unew[2:] = u[2:] - nu*(1.5*u[2:] - 2.0*u[1:-1] + 0.5*u[0:-2]) \
              +0.5*nu**2*(u[0:-2] - 2.0*u[1:-1] + u[2:])
  return unew

def solve(a, N, cfl, Tf):
    xmin, xmax = 0.0, 2*np.pi

    h = (xmax - xmin)/N
    dt= cfl * h / np.abs(a)
    nu= a * dt / h

    x = np.linspace(xmin, xmax, N+1)
    u = uinit(x)
    v = uinit(x)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    line1, = ax.plot(x, u, 'r') #Backward lw
    line2, = ax.plot(x, u, 'b') #Central lw
    line3, = ax.plot(x, u, 'y') #Exact Solution
    ax.set_xlabel('x'); ax.set_ylabel('u')
    plt.grid(True)
    plt.legend(('Backward Lax-Wendroff','Central Lax-Wendroff', 'Exact Solution'))
    plt.title('N='+str(N)+', CFL='+str(cfl))
    plt.draw(); plt.pause(0.1)
    wait = input("Press enter to continue ")

    t, it = 0.0, 0
    while t < Tf:
        u = update_back_lw(nu,u)
        v = update_lw(nu, v)
        t += dt; it += 1
        line1.set_ydata(u)
        line2.set_ydata(v)
        line3.set_ydata(uinit(x-a*t))
        plt.draw(); plt.pause(0.1)
    plt.show()

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of cells', default=100)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.98)
parser.add_argument('-a', type=float, help='Advection speed', default=1.0)
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)
parser.add_argument('-k', type = float, help='Frequency', default='1.0')
args = parser.parse_args()

# Run the solver
def uinit(x):
  return np.sin(args.k*x)
solve(args.a, args.N, args.cfl, args.Tf)