# Solver for Riemann problem of Linear Advection
# User chooses Ul, Ur

push!(LOAD_PATH,".")
using FV
using Grid
using EqLinAdv
using LinearAlgebra

grid_size = 50 # number of cells

# Plan - To not have to worry about periodicity, we can temporarily
# specify the boundary points to be the Dirichlet values.
# this will work until the wave moves to the boundary, so we will
# test till that short time.

xmin, xmax   = 0.0, 1.0 # domain
nvar         = 3        # number of variables
fprime(U,x,eq)  = [1.0 2.0 3.0;
                   0.0 -2.0 3.0;
                   0.0 0.0 3.0] # The matrix given by F'(U)
                # The last argument is absolutely dummy to solve
                # Euler with the same code. Should we make `eq` an
                # optional argument?
final_time = 1.0

num_flux   = upwind

# initial condition
# Ul, Ur       = [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]
# initial_value(x) = (x <= 0) ? Ul : Ur
boundary_value(x) = 0.0 # Dummy
boundary_condition = "Periodic"
initial_value(x) = [sin(2.0*pi*x),0.5*sin(2.0*pi*x),1.5*sin(2.0*pi*x)]

save_time_interval = 0.0

cfl = 0.0

# -----------------------------------
equation = get_equation(fprime)
problem = Problem((xmin,xmax), nvar, initial_value, boundary_value,
                  boundary_condition, final_time)
param = Parameters(grid_size, cfl, save_time_interval)
scheme = Scheme(num_flux)
p, anim = solve(equation, problem, scheme, param)
savefig(p, "final_soln.pdf")
gif(anim, "soln.gif") # would have been better in the solve function
                      # here because of VS Code

# TODO - compare with characteristic pictures in Ch 3
