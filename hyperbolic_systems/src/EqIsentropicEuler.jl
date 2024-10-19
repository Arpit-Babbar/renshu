module EqIsentropicEuler

using LinearAlgebra
using Plots
using DelimitedFiles
using StaticArrays
using LaTeXStrings
using UnPack

import ..FV1D: get_node_vars

# For efficiency, euler should also contain γm1,γm3,3γm1_2. That just might a
# kill of readability though. Also, we don't know whether accessing that
# far-away gamma is even better than computing it every time
struct IsentropicEuler
   γ::Float64
   eps::Float64
end

# Would have liked to solver

id_(x) = x

# pure function
function flux(x, U, eq::IsentropicEuler)
   gamma = eq.γ
   @unpack gamma = eps
   rho, rho_v = U        # density
   v = rho_v / rho       # velocity
   p = rho^gamma         # pressure
   F = SVector(rho, rho * v^2 + p / eps^2) # flux
   return F
end

# TODO - Find the best version by counting the number of operations!!

# The matrix fprime(U)
function fprime(x, U, eq::IsentropicEuler)
   return nothing # NOT USING THIS!
end

# function converting primitive variables to PDE variables
function primitive2pde(prim, γ) # primitive, viscosity
   rho, v = prim
   U = SVector(rho, rho * v)
   return U
end

# function converting pde variables to primitive variables
function pde2primitive(U, γ)
   rho, rho_v = U
   return SVector(rho, rho_v / rho)
end

function pressure(equation::IsentropicEuler, U)
    gamma = equation.γ
    p = U[1]^gamma
   return p
end

function get_node_vars(U, eq::IsentropicEuler, indices...)
    SVector(ntuple(v -> U[v, indices...], 2))
 end

#-------------------------------------------------------------------------------
# Numerical Fluxes
#-------------------------------------------------------------------------------
function lax_friedrich!(equation::IsentropicEuler, lam, Ul, Ur, x, Uf) # Numerical flux of face at x
   return nothing
end

function rusanov!(equation::IsentropicEuler, lam, Ul, Ur, x, Uf) # Numerical flux of face at x
   gamma  = equation.γ
   @unpack eps = equation
   rhol, ul = pde2primitive(Ul, gamma)
   rhor, ur = pde2primitive(Ur, gamma)
   pl, pr = pressure(equation, Ul), pressure(equation, Ur)

   cl, cr = sqrt(gamma*pl/rhol) / eps, sqrt(gamma*pr/rhor) / eps # sound speed
   λ = maximum(abs.((ul, ul-cl, ul+cl, ur, ur-cr, ur+cr))) # local wave speed
   Fl, Fr = flux(x, Ul, equation), flux(x, Ur, equation)
   F = 0.5*(Fl+Fr) - 0.5*λ*(Ur - Ul)
   return F
end

function steger_warming!(equation::IsentropicEuler, lam, Ul, Ur, x, Uf)
   return nothing
end

function roe!(equation::IsentropicEuler, lam, Ul, Ur, x, Uf)
   return nothing
end

function vanleer!(equation::IsentropicEuler, lam, Ul, Ur, x, Uf)
   return nothing
end

function hll!(equation::IsentropicEuler, lam, Ul, Ur, x, Uf)
   return nothing
end


function hllc!(equation::IsentropicEuler, lam, Ul, Ur, x, Uf)
   return nothing
end

numfluxes = Dict("lax_friedrich"  => lax_friedrich!,
                 "rusanov"        => rusanov!,
                 "steger_warming" => steger_warming!,
                 "vanleer"        => vanleer!,
                 "roe"            => roe!,
                 "hll"            => hll!,
                 "hllc"           => hllc!
                 )

#-------------------------------------------------------------------------------
# Plotting Functions
#-------------------------------------------------------------------------------
function initialize_plot(grid, problem, equation::IsentropicEuler, scheme, U)
   anim = Animation()
   xc = grid.xc
   nx = grid.nx
   @unpack nvar = problem
   Up = copy(U)
   for j=1:nx
      @views Up[:, j] = pde2primitive(U[:,j],equation.γ)
   end
   # Adding title as a subplot in itself
   p_title = plot(title = "$nx points, time = 0",
                          grid = false,
                          showaxis = false, bottom_margin = 0Plots.px)
   ymin, ymax = minimum(Up[1,1:nx]), maximum(Up[1,1:nx])
   p1 = @views plot(xc, Up[1,1:nx], legend=false, label = nothing,
                    ylim = (ymin-0.1, ymax+0.1), linestyle = :dot, color = :red,
                    markerstrokestyle = :dot, markershape = :circle,
                    markersize = 3, markerstrokealpha = 0)
   xlabel!(p1, "x"); ylabel!(p1, "Density")

   ymin, ymax = minimum(Up[2,1:nx]), maximum(Up[2,1:nx])
   p2 = @views plot(xc, Up[2,1:nx], legend=false, label = nothing,
                    ylim = (ymin-0.1, ymax+0.1), linestyle = :dot, color = :red,
                    markerstrokestyle = :dot, markershape = :circle,
                    markersize = 3, markerstrokealpha = 0)
   xlabel!(p2, "x"); ylabel!(p2, "Velocity")

   l = @layout[ a{0.01h}; b c]
   p = plot(p_title, p1, p2, layout = l) # Can this be nvar independent?
   frame(anim)
   return p, anim
end

function update_plot!(grid, problem, equation::IsentropicEuler, scheme, U, t, it, param, plt_data)
   p, anim = plt_data
   @unpack save_time_interval = param
   @unpack final_time = problem
   if save_time_interval > 0.0
      k1, k2 = ceil(t/save_time_interval), floor(t/save_time_interval)
      if !(abs(t-k1*save_time_interval) < 1e-10 ||
           abs(t-k2*save_time_interval) < 1e-10 ||
           abs(t-final_time) < 1e-10)
         return nothing
      end
   end
   xc = grid.xc
   nx = grid.nx
   @unpack nvar = problem
   Up = copy(U)
   for j=1:nx
      @views Up[:, j] = pde2primitive(U[:,j],equation.γ)
   end
   time = round(t, digits=3)
   title!(p[1], "$nx points, time = $time")
   for i=1:nvar
      y_lims = (minimum(Up[i,:])-0.1, maximum(Up[i,:])+0.1)
      ylims!(p[i+1],y_lims) # Bad solution
      p[i+1][1][:y] = @views Up[i,1:nx]
   end
   frame(anim)
end

function final_plot(plt_data, equation::IsentropicEuler)
   p, anim = plt_data
   soln_data = readdlm("toro_user_exact.dat", skipstart = 9);
   @views x = soln_data[:,1];
   @views dens_exact = soln_data[:,2];
   @views velx_exact = soln_data[:,3];
   plot!(p[2],x,dens_exact, label = nothing, color = :blue, legend=false)
   plot!(p[3],x,velx_exact, label = "Exact", color = :blue, legend=true)
   savefig("sol.png")
   return p
end

empty_func(x...)=nothing

function get_plot_funcs(skip_plotting)
   if skip_plotting == true
      return Dict("initialize_plot" => empty_func,
                  "update_plot!"    => empty_func,
                  "final_plot"      => empty_func)
   else
      return Dict("initialize_plot" => initialize_plot,
                  "update_plot!"    => update_plot!,
                  "final_plot"      => final_plot)
   end
end

# get_equation(γ) = Dict( "eq"              => IsentropicEuler(γ),
#                         "flux"            => flux,
#                         "fprime"          => fprime,
#                         "numfluxes"       => numfluxes,
#                         "name"            => "1D IsentropicEuler equations")
get_equation(γ, eps) = (; eq = IsentropicEuler(γ, eps), flux, fprime, numfluxes,
                     name = "1D Isentropic Euler equations")

export roe!
export lax_friedrich!
export hll!
export rusanov!
export steger_warming!
export vanleer!
export get_equation
export get_plot_funcs
export primitive2pde
export pde2primitive


end

# idx = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x