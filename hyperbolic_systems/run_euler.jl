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

grid_size = 50 # number of cells

# To not have to worry about periodicity, we can temporarily
# specify the boundary points to be the Dirichlet values.
# this will work until the wave moves to the boundary, so we will
# test till that short time.

xmin, xmax   = 0.0, 1.0 # domain
nvar         = 3        # number of variables
final_time   = 0.1
γ            = 1.4      # gas constant
disc_x = 0.5            # location of initial discontinuity

num_flux   = lax_friedrich # TODO -  Change to string

# initial condition
# Specify Ul, Ur in primitive coordinates
primitive_l, primitive_r  = [1.0, 0.0, 1.0], [0.125, 2.0, 1.1]
                            # density, velocity, pressure

Ul, Ur = primitive2pde(primitive_l, γ), primitive2pde(primitive_r, γ)
initial_value(x) = (x <= disc_x) ? Ul : Ur
boundary_value(x) = 0.0 # Dummy
boundary_condition = "Dirichlet"
# initial_value(x) = [sin(2.0*pi*x),sin(2.0*pi*x),sin(2.0*pi*x)]

save_time_interval = 0.0

cfl = 0.0

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
problem = Problem((xmin,xmax), nvar, initial_value, boundary_value,
                  boundary_condition, final_time)
param = Parameters(grid_size, cfl, save_time_interval)
scheme = Scheme(num_flux)
p, anim = solve(equation, problem, scheme, param)

#------------------------------------------------------------------------------
# Plot generated data
#------------------------------------------------------------------------------

soln_data = readdlm("toro_user_exact.dat", skipstart = 9);
@views x = soln_data[:,1];
@views dens_exact = soln_data[:,2];
@views pres_exact = soln_data[:,3];
@views velx_exact = soln_data[:,4];
plot!(p[2],x,dens_exact, label = nothing)
plot!(p[4],x,pres_exact, label = nothing)
plot!(p[3],x,velx_exact, label = "Exact", legend=true)
savefig(p, "final_soln.pdf")
gif(anim, "soln.gif") # would have been better in the solve function
                     # here because of VS Code
# TODO - compare with characteristic pictures in Ch 3

