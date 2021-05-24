module EqEuler

using LinearAlgebra
using PyPlot
using DelimitedFiles

struct Euler
   γ::Float64
end

# Would have liked to solver 

# pure function
function flux(x, U, eq::Euler) 
   ρ = U[1]        # density
   u = U[2] / U[1] # velocity
   E = U[3]        # energy
   p = (eq.γ - 1.0) * (E - 0.5 * ρ * u^2) # pressure
   F = [U[2], p + ρ * u^2, (E+p) * u] # flux
   return F
end
# TODO - Find the best version by counting the number of operations!!

# The matrix fprime(U)
function fprime(x, U, eq::Euler) 
   ρ = U[1]        # density
   u = U[2] / U[1] # velocity
   E = U[3]        # energy

   p = (eq.γ - 1.0) * (E - 0.5 * ρ * u^2) # pressure

   H = (E+p)/ρ

   A = [0.0                          1.0                  0.0;
        0.5*(eq.γ-3.0)*u^2       (3.0-eq.γ)*u     eq.γ-1.0;
        u*(0.5*(eq.γ-1.0)*u^2-H) H-(eq.γ-1.0)*u^2 eq.γ*u]
   return A
end

# different arguments from fprime in LinAdv. Can this be avoided?

function lax_friedrich(equation, lam, Ul, Ur, x) # Numerical flux of face at x
   eq = equation["eq"]
   Fl, Fr = flux(x, Ul, eq), flux(x, Ur, eq)
   value  = 0.5*(Fl+Fr) - 0.5 * lam * (Ur - Ul)
   return value
end

# function converting primitive variables to PDE variables
function primitive2pde(prim, γ) # primitive, viscosity
   U = [prim[1], prim[1]*prim[2], prim[3]/(γ-1.0) + prim[1]*prim[2]^2/2.0]
      # ρ    ,     ρ*u     ,        p/(γ-1.0) +     ρ*u^2/2.0
   return U
end

# function converting pde variables to primitive variables
function pde2primitive(U, γ)
   prim = [U[1], U[2]/U[1], (γ-1.0)*(U[3]-U[2]^2/(2.0*U[1]))]
   # prim=[ρ , u        , p]
   return prim
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
   plt.clf()
   xc = grid.xc
   nx = grid.nx
   eq = equation["eq"]
   Up = copy(U)
   for j=1:nx
      @views Up[:, j] = pde2primitive(U[:,j],eq.γ)
   end
   suptitle("Iteration $it, time $t")
   subplot(131)
   @views plot(xc, Up[1,1:nx])
   xlabel("x")
   ylabel("Density")
   subplot(132)
   @views plot(xc, Up[2,1:nx], "o", markersize=3)
   xlabel("x")
   ylabel("Velocity")
   subplot(133)
   @views plot(xc, Up[3,1:nx], "o", markersize=3)
   xlabel("x")
   ylabel("Pressure")
   plt.pause(0.1)
end

function plot_final_soln(grid, equation, problem, U, t, it, param)
   close("all")
   soln_data = readdlm("toro_user_exact.dat", skipstart = 9)
   @views x = soln_data[:,1]
   @views dens_exact = soln_data[:,2]
   @views pres_exact = soln_data[:,3]
   @views velx_exact = soln_data[:,4]
   xc = grid.xc
   nx = grid.nx
   eq = equation["eq"]
   Up = copy(U)
   for j=1:nx
      @views Up[:, j] = pde2primitive(U[:,j],eq.γ)
   end
   subplots(1,3)
   suptitle("Iteration $it, time $t")
   subplot(131)
   @views plot(xc, Up[1,1:nx], "o", markersize=3)
   plot(x, dens_exact)
   legend(("Numerical", "Exact"))
   xlabel("x")
   ylabel("Density")
   subplot(132)
   @views plot(xc, Up[2,1:nx], "o", markersize=3)
   plot(x, velx_exact)
   legend(("Numerical", "Exact"))
   xlabel("x")
   ylabel("Velocity")
   subplot(133)
   @views plot(xc, Up[3,1:nx], "o", markersize=3)
   plot(x, pres_exact)
   legend(("Numerical", "Exact"))
   xlabel("x")
   ylabel("Pressure")
end
# TODO - Make a plot function that does the common parts in the function. Something like
# fig, ax = Init_triple_plot()
# add_to_figure(ax) <- This looks like something that won't freaking work!!
get_equation(γ) = Dict( "eq"              => Euler(γ),
                        "flux"            => flux,
                        "fprime"          => fprime,
                        "plot_solution"   => plot_solution,
                        "plot_final_soln" => plot_final_soln,
                        "name"            => "1D Euler equations")


export lax_friedrich
export get_equation
export primitive2pde
export pde2primitive
export plot_final_soln

end
