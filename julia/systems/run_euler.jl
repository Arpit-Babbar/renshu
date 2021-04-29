# Solver for Riemann problem of Linear Advection
# User chooses Ul, Ur

push!(LOAD_PATH,".")
using FV
using Grid
using EqEuler
using LinearAlgebra

grid_size = 300 # number of cells

# Plan - To not have to worry about periodicity, we can temporarily
# specify the boundary points to be the Dirichlet values.
# this will work until the wave moves to the boundary, so we will
# test till that short time.

xmin, xmax   = 0.0, 1.0 # domain
nvar         = 3        # number of variables
final_time   = 0.20
gamma  = 1.4
disc_x = 0.5

num_flux   = lax_friedrich

# initial condition
# Specify Ul, Ur in primitive coordinates and make a function that
# converts primitives to the usual format and another function that
# brings them back. The other function will be used for plotting
primitive_l, primitive_r  = [1.0, 0.0, 1.0], [0.125, 0.0, 0.1]
                            # density, velocity, pressure
println("Solving Riemann problem for Euler equation with following parameters -")
println("Discontinuity of initial data is at", disc_x)
println("gamma = ", gamma)
println("Tf    = ", final_time)
println("Initial Condition L/R of [dens, vel, pres] is")
show(stdout, primitive_l)
println("")
show(stdout, primitive_r)
println("")

println("First we solve it exactly")
p_l_s, p_r_s = array2string(primitive_l), array2string(primitive_r)
run(`python ./ToroExact/toro_exact.py -p user -l $p_l_s -r $p_r_s -x $disc_x -t $final_time`)

Ul, Ur = primitive2pde(primitive_l, gamma), primitive2pde(primitive_r, gamma)

initial_value(x) = (x <= disc_x) ? Ul : Ur
boundary_value(x) = 0.0 # Dummy
boundary_condition = "Dirichlet"
# initial_value(x) = [sin(2.0*pi*x),sin(2.0*pi*x),sin(2.0*pi*x)]

save_time_interval = 0.0

cfl = 0.0

# -----------------------------------
equation = get_equation(gamma)
problem = Problem((xmin,xmax), nvar, initial_value, boundary_value,
                  boundary_condition, final_time)
param = Parameters(grid_size, cfl, save_time_interval)
scheme = Scheme(num_flux)
solve(equation, problem, scheme, param)

# TODO - compare with characteristic pictures in Ch 3
