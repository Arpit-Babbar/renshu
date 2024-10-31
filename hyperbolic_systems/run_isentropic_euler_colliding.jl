# Solver for Riemann problem of Linear Advection
# User chooses Ul, Ur

using FV1D
using FV1D.FV
using FV1D.Grid
using FV1D.EqIsentropicEuler
using LinearAlgebra
using StaticArrays
using DelimitedFiles
using Plots
using UnPack
gr(size = (750, 565)) # use plotly for interactivity

grid_size = 80 # number of cells

# To not have to worry about periodicity, we can temporarily
# specify the boundary points to be the Dirichlet values.
# this will work until the wave moves to the boundary, so we will
# test till that short time.

xmin, xmax   = -1.0, 1.0 # domain
nvar         = 2        # number of variables
final_time   = 0.04
γ            = 1.4      # gas constant

numflux = "flux_central"

function initial_value!(U, x, equation)
   @unpack epsilon = equation
   @show epsilon
   gamma = equation.γ
   rho = 0.955 + 0.5*epsilon*(1.0-cospi(2.0*x))
   vel = -sign(x) * sqrt(gamma) * (1.0 - cospi(2.0*x))
   U .= rho, rho * vel
   return nothing
end

# Also the exact solution
function boundary_value(x, t, equation)
   @unpack epsilon = equation
   gamma = equation.γ
   rho = 0.955 + 0.5*epsilon*(1.0-cospi(2.0*x))
   vel = -sign(x) * sqrt(gamma) * (1.0 - cospi(2.0*x))
   return SVector(rho, rho * vel)
end
boundary_condition = "Periodic"

save_time_interval = 0.01*final_time
skip_plotting      = false
plotters = EqIsentropicEuler.get_plot_funcs(skip_plotting)

cfl = 0.0
Ccfl = 0.9

#------------------------------------------------------------------------------
# FVM Solver
#------------------------------------------------------------------------------
epsilon = 0.1
equation = EqIsentropicEuler.get_equation(γ, epsilon)
initial_value!(U, x) = initial_value!(U, x, equation.eq)
boundary_value(x, t) = boundary_value(x, t, equation.eq)
problem = Problem((xmin,xmax), nvar, initial_value!, boundary_value,
                  boundary_condition, final_time)
param = Parameters(grid_size, cfl, Ccfl, save_time_interval)
scheme = Scheme(equation, numflux)
@time sol = solve(equation, problem, scheme, param, x->x, plotters)

@show sol.l2
@show sol.l1
@show sol.linf

@unpack p, anim = sol
savefig(p, "final_soln.png")
gif(anim, "soln.mp4", fps = 1) # would have been better in the solve function
                     # here because of VS Code
