module EqLinAdv

using LinearAlgebra
using PyPlot

struct LinAdv
   fprime::Function
end

flux(x, U, eq::LinAdv) = eq.fprime(U,x) * U
                        # matrix multiplication
# Doubt - lax_friedrich shouldn't work at the right Dirichlet bc, right?
# Over there, we can't apply the Dirichlet condition.
# Yes, that's a bug. Must fix!!
function lax_friedrich(equation, lam, Ul, Ur, x) # Numerical flux of face at x
   eq = equation["eq"]
   Fl, Fr = flux(x, Ul, eq), flux(x, Ur, eq)
   value  = 0.5*(Fl+Fr) - 0.5 * lam * (Ur - Ul)
   return value
end

# Looks expensive for a non-constant speed. Can it be fixed?
function upwind(equation, lam, Ul, Ur, x) # Numerical flux of face at x
   # eigen_decomp   = eigen(eq.fprime(0.5 * (Ul+Ur), x, eq)) # Roe's scheme had (f(Ul)-f(Ur))/(Ul-Ur)
   #                               # appearing hard to extend to high dimensions because Ul, Ur
   #                               # are now vectors and can't be put in denominators.
   #                               # for linear advection, it is working because fprime is a constant matrix.
   eigen_decomp = equation["eigen_decomp"]
   evalues, evecs = eigen_decomp.values, eigen_decomp.vectors
   lp, lm     = max.(evalues, 0.0), min.(evalues, 0.0)
   lamp, lamm = diagm(lp), diagm(lm)
   Fp = evecs * lamp * inv(evecs) * Ul
   Fm = evecs * lamm * inv(evecs) * Ur
   value  = Fp + Fm
   return value
end

# Repetetive calculation of eigen_decomp every time 
# you compute the exact solution. Should try to fix
# The issue with fixing is that
function compute_exact_soln!(grid, equation, problem, t, Ue)
   nx = grid.nx
   xc = grid.xc
   initial_value = problem["initial_value"]
   eigen_decomp = equation["eigen_decomp"]
   nvar = problem["nvar"]
   Lam, eigen_vecs = eigen_decomp.values, eigen_decomp.vectors
   for j=1:nx
      for i=1:nvar
         # Room for simplification? Yes, too many brackets
         Ue[i,j] = (inv(eigen_vecs) * initial_value(xc[j] - Lam[i] * t))[i]
      end
   end
   for j=1:nx
      Ue[:,j] = eigen_vecs * Ue[:,j]
   end
   return nothing
end

# Fix - Make the plot_solution a part of PDE. And, in the plot function
# compute exact solution. And, create a boundary_value function for Dirichlet bc
# function plot_solution(grid, U, Ue, t, it, param)
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
   nvar = problem["nvar"]
   Ue = zeros(nvar, nx)# Need to move it elsewhere to avoid
                # computing it every time
   compute_exact_soln!(grid, equation, problem, t, Ue)
   suptitle("Iteration $it, time $t")
   subplot(131)
   @views plot(xc, U[1,1:nx])
   @views plot(xc, Ue[1,:])
   legend(("Approximate", "Exact"))
   xlabel("x")
   ylabel("\$U_1\$")
   subplot(132)
   @views plot(xc, U[2,1:nx])
   @views plot(xc, Ue[2,:])
   legend(("Approximate", "Exact"))
   xlabel("x")
   ylabel("\$U_2\$")
   subplot(133)
   @views plot(xc, U[3,1:nx])
   @views plot(xc, Ue[3,:])
   legend(("Approximate", "Exact"))
   xlabel("x")
   ylabel("\$U_3\$")
   plt.pause(0.1)
end


get_equation(fprime) = Dict("eq"                  => LinAdv(fprime),
                            "flux"                => flux,
                            "fprime"              => fprime,
                            "eigen_decomp"        => eigen(fprime(0.0,0.0,0.0)), # dummy inputs
                            # Since we are solving constant variables advection, it is best to compute
                            # eigen_decomp only once. The moment we go to variables advection, this would Need
                            # to be removed
                            "compute_exact_soln!" => compute_exact_soln!,
                            "plot_solution"       => plot_solution,
                            "plot_final_soln"     => plot_solution,
                            "name"                => "Linear advection equation")

export LinAdv # To define fprime in run file
export lax_friedrich
export upwind
export get_equation
export plot_solution

end