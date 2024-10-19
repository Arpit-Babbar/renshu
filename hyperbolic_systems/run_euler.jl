# Solver for Riemann problem of Linear Advection
# User chooses Ul, Ur

using FV1D
using FV1D.FV
using FV1D.Grid
using FV1D.EqEuler
using LinearAlgebra
using DelimitedFiles
using StaticArrays
using Plots
using UnPack
gr(size = (750, 565)) # use plotly for interactivity

grid_size = 1000 # number of cells

# To not have to worry about periodicity, we can temporarily
# specify the boundary points to be the Dirichlet values.
# this will work until the wave moves to the boundary, so we will
# test till that short time.

xmin, xmax   = 0.0, 1.0 # domain
nvar         = 3        # number of variables
final_time   = 0.2
γ            = 1.4      # gas constant
disc_x = 0.3            # location of initial discontinuity

numflux = "hll"

id(x) = x

# Choose between id, CuArray
promoter = x -> x

# initial condition
# Specify Ul, Ur in primitive coordinates
                            # density, velocity, pressure

function prim2con(prim, γ) # primitive, gas constant
   U = SVector{3,Float64}(prim[1], prim[1]*prim[2], prim[3]/(γ-1.0) + 0.5*prim[1]*prim[2]^2)
   #           ρ    ,     ρ*u     ,        p/(γ-1.0) +     ρ*u^2/2.0
   return U
end

function prim2con!(u, v, γ)
   u .= prim2con(v, γ)
   return nothing
end

function initial_value_general!(x, disc_x)
   primitive_l, primitive_r  = (1.0, 0.75, 1.0), (0.125, 0.0, 0.1)

   γ            = 1.4      # gas constant
   if x <= disc_x
      return prim2con(primitive_l, γ)
   else
      return prim2con(primitive_r, γ)
   end
   return nothing
end

function initial_value_general!(U, x, disc_x)
   U .= initial_value_general!(x, disc_x)
   return nothing
end
# TODO - This gives dynamic invocation error.
initial_value!(U, x) = initial_value_general!(U, x, disc_x)

boundary_value(x, t) = 0.0 # Dummy
boundary_condition = "Dirichlet"
# initial_value(x) = [sin(2.0*pi*x),sin(2.0*pi*x),sin(2.0*pi*x)]

save_time_interval = 0.1*final_time
skip_plotting      = false
plotters = get_plot_funcs(skip_plotting)

cfl = 0.0
Ccfl = 0.9

#------------------------------------------------------------------------------
# Print parameters to screen
#------------------------------------------------------------------------------
primitive_l = pde2primitive(initial_value_general!(-Inf, disc_x), γ)
primitive_r = pde2primitive(initial_value_general!(Inf, disc_x), γ)
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
# print("python3 ./ToroExact/toro_exact.py -p user -l $p_l_s -r $p_r_s -x $disc_x -t $final_time")

#------------------------------------------------------------------------------
# FVM Solver
#------------------------------------------------------------------------------
equation = get_equation(γ)
problem = Problem((xmin,xmax), nvar, initial_value!, boundary_value,
                  boundary_condition, final_time)
param = Parameters(grid_size, cfl, Ccfl, save_time_interval)
scheme = Scheme(equation, numflux)
@time sol = solve(equation, problem, scheme, param, promoter, plotters)

@unpack p, anim = sol
savefig(p, "final_soln.png")
gif(anim, "soln.gif", fps = 1) # would have been better in the solve function
                     # here because of VS Code

