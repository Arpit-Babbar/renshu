# Solver for Riemann problem of Linear Advection
# User chooses Ul, Ur

using FV1D
using FV1D.FV
using FV1D.Grid
using FV1D.EqEuler
using LinearAlgebra
using StaticArrays
using DelimitedFiles
using Plots
using UnPack
gr(size = (750, 565)) # use plotly for interactivity

grid_size = 400 # number of cells

# To not have to worry about periodicity, we can temporarily
# specify the boundary points to be the Dirichlet values.
# this will work until the wave moves to the boundary, so we will
# test till that short time.

xmin, xmax   = 0.0, 1.0 # domain
nvar         = 3        # number of variables
final_time   = 1.0
γ            = 1.4      # gas constant
disc_x = 0.3            # location of initial discontinuity

numflux = "rusanov"

id(x) = x

# Choose between id, CuArray
promoter = x -> x

# initial condition
# Specify Ul, Ur in primitive coordinates
primitive_l, primitive_r  = [1.0, 0.75, 1.0], [0.125, 0.0, 0.1]
                            # density, velocity, pressure

function prim2con(prim, γ) # primitive, gas constant
   U = SVector{3,Float64}(prim[1], prim[1]*prim[2], prim[3]/(γ-1.0) + 0.5*prim[1]*prim[2]^2)
   #           ρ    ,     ρ*u     ,        p/(γ-1.0) +     ρ*u^2/2.0
   return U
end

function prim2con!(u, v, γ)
   u[1] = v[1]
   u[2] = v[1]*v[2]
   u[3]  = v[3]/(γ-1.0) + 0.5*v[1]*v[2]^2
   return nothing
end

function initial_value_general!(U, x, disc_x)

   γ            = 1.4      # gas constant
   if x <= disc_x
      prim2con!(U, (1.0, 0.75, 1.0), γ)
   else
      prim2con!(U, (0.125, 0.0, 0.1), γ)
   end
   return nothing
end

function initial_value!(U, x)
   γ = 1.4
   rho = 1.0 + 0.5 * sinpi(2.0 * x)
   vel = 1.0
   p = 1.0
   U .= rho, rho * vel, p / (γ - 1.0) + 0.5 * rho * vel^2
   return nothing
end

# TODO - This gives dynamic invocation error.
# initial_value!(U, x) = initial_value_general!(U, x, disc_x)

function boundary_value(x, t)
   γ = 1.4
   rho = 1.0 + 0.5 * sinpi(2.0 * (x - t))
   vel = 1.0
   p = 1.0
   return SVector(rho, rho * vel, p / (γ - 1.0) + 0.5 * rho * vel^2)
   return nothing
end
boundary_condition = "Periodic"
# initial_value(x) = [sin(2.0*pi*x),sin(2.0*pi*x),sin(2.0*pi*x)]

save_time_interval = 0.1*final_time
skip_plotting      = false
plotters = get_plot_funcs(skip_plotting)

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

@show sol.l2
@show sol.l1
@show sol.linf

@unpack p, anim = sol
savefig(p, "final_soln.png")
gif(anim, "soln.mp4", fps = 1) # would have been better in the solve function
                     # here because of VS Code
