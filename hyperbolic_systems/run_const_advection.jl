# Solver for Riemann problem of Linear Advection
# User chooses Ul, Ur

using FV1D
using FV1D.FV
using FV1D.Grid
using FV1D.EqLinAdv
using StaticArrays
using LinearAlgebra
using Plots

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
fprime(U,x,eq)  = [1.0 0.0 0.0;
                   0.0 1.0 0.0;
                   0.0 0.0 1.0] # The matrix given by F'(U)
                # The last argument is absolutely dummy to solve
                # Euler with the same code. Should we make `eq` an
                # optional argument?
final_time = 1.0

numflux   = "upwind"

# initial condition
# Ul, Ur       = [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]
# initial_value(x) = (x <= 0) ? Ul : Ur
boundary_value(x, t) = initial_value(x) # Dummy

function boundary_value(x,t,problem,equation)
   eigen_decomp = equation.eigen_decomp
   nvar = problem.nvar
   Lam, eigen_vecs = eigen_decomp.values, eigen_decomp.vectors
   Ue = zeros(3)
   for i=1:3
      # Room for simplification? Yes, too many brackets
      Ue[i] = (inv(eigen_vecs) * initial_value(x - Lam[i] * t))[i]
   end
   Ue .= eigen_vecs * Ue
   return Ue
end

boundary_value(x,t) = boundary_value(x,t,problem, equation)

boundary_condition = "Periodic"
function initial_value(x)
   return SVector(sin(2.0*pi*x),0.5*sin(2.0*pi*x),1.5*sin(2.0*pi*x))
end
function initial_value(U, x)
   U .= SVector(x)
end

save_time_interval = final_time
skip_plotting = false
cfl = 0.0 # Currently not used
Ccfl = 0.9

# -----------------------------------
plotters = EqLinAdv.get_plot_funcs(skip_plotting)
equation = EqLinAdv.get_equation(fprime)
problem = Problem((xmin,xmax), nvar, initial_value, boundary_value,
                  boundary_condition, final_time)
param = Parameters(grid_size, cfl, Ccfl, save_time_interval)
scheme = Scheme(equation, numflux)
@time sol = solve(equation, problem, scheme, param, u -> u, plotters)


@unpack p, anim = sol
savefig(p, "final_soln.png")
gif(anim, "soln.gif", fps = 1) # would have been better in the solve function
                     # here because of VS Code
plot(p, legend=true) # final solution
