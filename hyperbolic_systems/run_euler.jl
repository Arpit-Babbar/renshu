# Solver for Riemann problem of Linear Advection
# User chooses Ul, Ur

push!(LOAD_PATH,".")
using FV
using Grid
using EqEuler
using LinearAlgebra
using DelimitedFiles
using Plots
plotly(size = (750, 565)) # use plotly for interactivity

grid_size = 100 # number of cells

# To not have to worry about periodicity, we can temporarily
# specify the boundary points to be the Dirichlet values.
# this will work until the wave moves to the boundary, so we will
# test till that short time.

xmin, xmax   = 0.0, 1.0 # domain
nvar         = 3        # number of variables
final_time   = 0.2
γ            = 1.4      # gas constant
disc_x = 0.3            # location of initial discontinuity

numflux   = "hll"

# initial condition
# Specify Ul, Ur in primitive coordinates
primitive_l, primitive_r  = [1.0, 0.75, 1.0], [0.125, 0.0, 0.1]
                            # density, velocity, pressure

Ul, Ur = primitive2pde(primitive_l, γ), primitive2pde(primitive_r, γ)
initial_value(x) = (x <= disc_x) ? Ul : Ur
boundary_value(x) = 0.0 # Dummy
boundary_condition = "Dirichlet"
# initial_value(x) = [sin(2.0*pi*x),sin(2.0*pi*x),sin(2.0*pi*x)]

save_time_interval = final_time
skip_plotting      = false

cfl = 0.0
Ccfl = 0.9

#------------------------------------------------------------------------------
# Print parameters to screen
#------------------------------------------------------------------------------
println("Solving Riemann problem for Euler equation with following parameters -")
println("Discontinuity of initial data is at", disc_x)
println("gamma = ", γ)
println("Tf    = ", final_time)
println("Initial Condition L/R of [dens, vel, pres] is")
show(stdout, primitive_l)
println("")
show(stdout, primitive_r)
println("")

#------------------------------------------------------------------------------
# Exact Solver
#------------------------------------------------------------------------------
println("Solving exactly for final time")
p_l_s, p_r_s = array2string(primitive_l), array2string(primitive_r)
run(`python3 ./ToroExact/toro_exact.py -p user -l $p_l_s -r $p_r_s -x $disc_x -t $final_time`)

#------------------------------------------------------------------------------
# FVM Solver
#------------------------------------------------------------------------------
equation = get_equation(γ)
plotters = get_plot_funcs(skip_plotting)
problem = Problem((xmin,xmax), nvar, initial_value, boundary_value,
                  boundary_condition, final_time)
param = Parameters(grid_size, cfl, Ccfl, save_time_interval)
scheme = Scheme(equation, numflux)
@time plt_data = solve(equation, problem, scheme, param, plotters)
#------------------------------------------------------------------------------
# Plot generated data
#------------------------------------------------------------------------------
if plt_data !== nothing
   p, anim = plt_data
   gif(anim, "soln.gif", fps = 5) # would have been better in the solve function
   # here because of VS Code
   plot(p, legend=true) # final solution
end
# TODO - compare with characteristic pictures in Ch 3

