module EqEuler

using LinearAlgebra
using PyPlot

struct Euler
   gamma::Float64
end

# Would have liked to solver 

# pure function
function flux(x, U, eq::Euler) 
   rho = U[1]        # density
   u   = U[2] / U[1] # velocity
   E   = U[3]        # energy
   p   = (eq.gamma - 1.0) * (E - 0.5 * rho * u^2) # pressure
   F   = [U[2], p + rho * u^2, (E+p) * u] # flux
   return F
end
# TODO - Find the best version by counting the number of operations!!

# The matrix fprime(U)
function fprime(x, U, eq::Euler) 
   rho = U[1]        # density
   u   = U[2] / U[1] # velocity
   E   = U[3]        # energy

   p   = (eq.gamma - 1.0) * (E - 0.5 * rho * u^2) # pressure

   H = (E+p)/rho

   A = [0.0                          1.0                  0.0;
        0.5*(eq.gamma-3.0)*u^2       (3.0-eq.gamma)*u     eq.gamma-1.0;
        u*(0.5*(eq.gamma-1.0)*u^2-H) H-(eq.gamma-1.0)*u^2 eq.gamma*u]
   return A
end

# different arguments from fprime in LinAdv. Can this be avoided?

function lax_friedrich(equation, lam, Ul, Ur, x) # Numerical flux of face at x
   eq = equation["eq"]
   Fl, Fr = flux(x, Ul, eq), flux(x, Ur, eq)
   value  = 0.5*(Fl+Fr) - 0.5 * lam * (Ur - Ul)
   return value
end

# Fix - Make the plot_solution a part of PDE. And, in the plot function
# compute exact solution. And, create a boundary_value function for Dirichlet bc
# function plot_solution(grid, U, Ue, t, it, param)
# Wouldn't it be better to make a general plot function in FV.jl and learn to add exact
# solution to that plot? Yeah, that'd be better.
function plot_solution(grid, equation, problem, U, t, it, param)
   save_time_interval = param["save_time_interval"]
   if save_time_interval > 0.0
      k1, k2 = ceil(t/save_time_interval), floor(t/save_time_interval)
      if (abs(t-k1*save_time_interval) < 1e-10 ||
          abs(t-k2*save_time_interval) < 1e-10)
         nothing
      else
         return nothing
      end
   end
   plt. clf()
   xc = grid.xc
   nx = grid.nx
   suptitle("Iteration $it, time $t")
   subplot(131)
   @views plot(xc, U[1,1:nx])
   legend(("Approximate", "Exact"))
   xlabel("x")
   ylabel("\$U_1\$")
   subplot(132)
   @views plot(xc, U[2,1:nx])
   legend(("Approximate", "Exact"))
   xlabel("x")
   ylabel("\$U_2\$")
   subplot(133)
   @views plot(xc, U[3,1:nx])
   legend(("Approximate", "Exact"))
   xlabel("x")
   ylabel("\$U_3\$")
   plt.pause(0.1)
end

get_equation(gamma) = Dict( "eq"            => Euler(gamma),
                            "flux"          => flux,
                            "fprime"        => fprime,
                            "plot_solution" => plot_solution,
                            "name"          => "1D Euler equations")

export lax_friedrich
export get_equation

end