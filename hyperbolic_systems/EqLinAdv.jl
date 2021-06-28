module EqLinAdv

using LinearAlgebra
# using PyPlot
using Plots
using LaTeXStrings

struct LinAdv
   fprime::Function
end

flux(x, U, eq::LinAdv) = eq.fprime(U, x, eq) * U
                        # matrix multiplication
# Doubt - lax_friedrich shouldn't work at the right Dirichlet bc, right?
# Over there, we can't apply the Dirichlet condition.
# Yes, that's a bug. Must fix!!
function lax_friedrich(equation, lam, Ul, Ur, x) # Numerical flux of face at x
   eq = equation["eq"]
   Fl, Fr = flux(x, Ul, eq), flux(x, Ur, eq)
   value  = 0.5*(Fl+Fr) - 0.5*lam*(Ur-Ul)
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

function initialize_plot(grid, problem, equation, scheme, U)
   anim = Animation()
   xc = grid.xc
   nx = grid.nx
   nvar = problem["nvar"]
   numflux = scheme["numflux_ind"]      # numflux as string
   numflux = replace(numflux, "_"=>" ") # Remove underscore
   numflux = titlecase(numflux)         # Capitalize
   # Adding title as a subplot in itself
   p_title = plot(title = "$numflux flux, $nx points, time = 0", grid = false,
                          showaxis = false, bottom_margin = 0Plots.px)
   p = [p_title]
   for i=1:nvar
      ymin, ymax = minimum(U[i,1:nx]), maximum(U[i,1:nx])
      p0 = @views plot(xc, U[i,1:nx], label="Approximate", ylim = (ymin-0.1, ymax+0.1))
      @views plot!(p0, xc, U[i,1:nx], label="Exact",ylim = (-0.1,1.1))
      xlabel!(p0, "x"); ylabel!(p0, "U")
      push!(p, p0)
   end
   l = @layout[ a{0.01h}; b c d]
   p = plot(p[1], p[2], p[3], p[4], layout = l) # Can this be nvar independent?
   frame(anim)
   return p, anim
end

function update_plot!(grid, problem, equation, scheme, U, t, it, param, p, anim)
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
   xc = grid.xc
   nx = grid.nx
   nvar = problem["nvar"]
   numflux = scheme["numflux_ind"]      # numflux as string
   numflux = replace(numflux, "_"=>" ") # Remove underscore
   numflux = titlecase(numflux)         # Capitalize
   Ue = zeros(nvar, nx)# Need to move it elsewhere to avoid
                # computing it every time
   compute_exact_soln!(grid, equation, problem, t, Ue)
   title!(p[1], "$numflux flux, $nx points, time = $time")
   for i=1:nvar
      y_lims = (min(minimum(Ue[i,:]),minimum(U[i,:]))-0.1,
                max(maximum(Ue[i,:]),maximum(U[i,:]))+0.1)
      ylims!(p[i+1],y_lims) # Bad approach
      p[i+1][1][:y] = @views U[i,1:nx]
   end
   frame(anim)
end

numfluxes = Dict("upwind"        => upwind,
                 "lax_friedrich" => lax_friedrich)

get_equation(fprime) = Dict("eq"                  => LinAdv(fprime),
                            "flux"                => flux,
                            "fprime"              => fprime,
                            "eigen_decomp"        => eigen(fprime(0.0,0.0,0.0)), # dummy inputs
                            # Since we are solving constant variables advection, it is best to compute
                            # eigen_decomp only once. The moment we go to variables advection, this would Need
                            # to be removed
                            "compute_exact_soln!" => compute_exact_soln!,
                            "initialize_plot"     => initialize_plot,
                            "update_plot!"        => update_plot!,
                            "numfluxes"           => numfluxes,
                            "name"                => "Linear advection equation")



export LinAdv # To define fprime in run file
export lax_friedrich
export upwind
export get_equation
export plot_solution

end