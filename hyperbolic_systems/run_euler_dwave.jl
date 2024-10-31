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

grid_size = 160 # number of cells

# To not have to worry about periodicity, we can temporarily
# specify the boundary points to be the Dirichlet values.
# this will work until the wave moves to the boundary, so we will
# test till that short time.

xmin, xmax   = 0.0, 1.0 # domain
nvar         = 3        # number of variables
final_time   = 0.1
γ            = 1.4      # gas constant
disc_x = 0.3            # location of initial discontinuity

numflux = "rusanov"

id(x) = x

# Choose between id, CuArray
promoter = x -> x

function initial_value!(U, x)
   γ = 1.4
   rho = 1.0 + 0.5 * sinpi(2.0 * x)
   vel = 1.0
   p = 1.0
   U .= rho, rho * vel, p / (γ - 1.0) + 0.5 * rho * vel^2
   return nothing
end

# Also the exact solution
function boundary_value(x, t)
   γ = 1.4
   rho = 1.0 + 0.5 * sinpi(2.0 * (x - t))
   vel = 1.0
   p = 1.0
   return SVector(rho, rho * vel, p / (γ - 1.0) + 0.5 * rho * vel^2)
end
boundary_condition = "Periodic"

save_time_interval = 0.1*final_time
skip_plotting      = false
plotters = EqEuler.get_plot_funcs(skip_plotting)

cfl = 0.0
Ccfl = 0.5

#------------------------------------------------------------------------------
# FVM Solver
#------------------------------------------------------------------------------
equation = EqEuler.get_equation(γ)
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
