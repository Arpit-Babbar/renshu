module EqIsentropicEuler

using LinearAlgebra
using Plots
using DelimitedFiles
using StaticArrays
using LaTeXStrings
using UnPack
using DelimitedFiles

import ..FV1D.FV: get_node_vars, compute_dt

import ..FV1D.EqEuler: initialize_plot, update_plot!, final_plot

# For efficiency, euler should also contain γm1,γm3,3γm1_2. That just might a
# kill of readability though. Also, we don't know whether accessing that
# far-away gamma is even better than computing it every time
struct IsentropicEuler
   γ::Float64
   epsilon::Float64
end

# Would have liked to solver

id_(x) = x

# pure function
function flux(x, U, eq::IsentropicEuler)
   gamma = eq.γ
   @unpack epsilon = eq
   rho, rho_v = U        # density
   v = rho_v / rho       # velocity
   p = rho^gamma         # pressure
   F = SVector(rho * v, rho * v^2 + p / epsilon^2) # flux
   return F
end

# TODO - Find the best version by counting the number of operations!!

# The matrix fprime(U)
function fprime(x, U, equation::IsentropicEuler)
   gamma = equation.γ
   @unpack epsilon = equation
   rho, u = pde2primitive(U, gamma)
   p = pressure(equation, U)

   c = sqrt(gamma*p/rho) / epsilon # sound speed
   λ = abs(u) + c # local wave speed
   return λ # NOT USING THIS!
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

 function compute_dt(equation::IsentropicEuler, fprime, scheme, param, grid, Ua, dt, lam)
   @unpack Ccfl = param
   nx, dx = grid.nx, grid.dx
   xc     = grid.xc
   lam_loc = 0.0
   dt_loc  = 1.0
   gamma = equation.γ
   @unpack epsilon = equation
   for i=1:nx
      ua   = Ua[:, i]
      rho, u = pde2primitive(ua, gamma)
      p = pressure(equation, ua)

      c = sqrt(gamma*p/rho) / epsilon # sound speed
      λ = abs(u) + c # local wave speed
      lam0 = λ
      lam_loc  = max(lam_loc, lam0)
      dt_loc   = min(dt_loc, dx[i]/lam0)
   end
   dt[1]   = dt_loc
   lam[1] = lam_loc[1]
   dt .= Ccfl*dt
end

#-------------------------------------------------------------------------------
# Numerical Fluxes
#-------------------------------------------------------------------------------
function lax_friedrich!(equation::IsentropicEuler, lam, Ul, Ur, x, Uf) # Numerical flux of face at x
   return nothing
end

function rusanov!(equation::IsentropicEuler, lam, Ul, Ur, x, Uf) # Numerical flux of face at x
   gamma  = equation.γ
   @unpack epsilon = equation
   rhol, ul = pde2primitive(Ul, gamma)
   rhor, ur = pde2primitive(Ur, gamma)
   pl, pr = pressure(equation, Ul), pressure(equation, Ur)

   cl, cr = sqrt(gamma*pl/rhol) / epsilon, sqrt(gamma*pr/rhor) / epsilon # sound speed
   λ = maximum(abs.((ul, ul-cl, ul+cl, ur, ur-cr, ur+cr))) # local wave speed
   Fl, Fr = flux(x, Ul, equation), flux(x, Ur, equation)
   F = 0.5*(Fl+Fr) - 0.5*λ*(Ur - Ul)
   return F
end

function flux_central(equation::IsentropicEuler, lam, Ul, Ur, x, Uf) # Numerical flux of face at x
   Fl, Fr = flux(x, Ul, equation), flux(x, Ur, equation)
   F = 0.5*(Fl+Fr)
   return F
end

@inline function ln_mean(x, y)
   epsilon_f2 = 1.0e-4
   f2 = (x * (x - 2 * y) + y * y) / (x * (x + 2 * y) + y * y) # f2 = f^2
   if f2 < epsilon_f2
       return (x + y) / @evalpoly(f2, 2, 2/3, 2/5, 2/7)
   else
       return (y - x) / log(y / x)
   end
end

@inline function stolarsky_mean(x, y, gamma)
   epsilon_f2 = 1.0e-4
   f2 = (x * (x - 2 * y) + y * y) / (x * (x + 2 * y) + y * y) # f2 = f^2
   if f2 < epsilon_f2
       # convenience coefficients
       c1 = (1 / 3) * (gamma - 2)
       c2 = -(1 / 15) * (gamma + 1) * (gamma - 3) * c1
       c3 = -(1 / 21) * (2 * gamma * (gamma - 2) - 9) * c2
       return 0.5 * (x + y) * @evalpoly(f2, 1, c1, c2, c3)
   else
       return (gamma - 1) / gamma * (y^gamma - x^gamma) /
              (y^(gamma - 1) - x^(gamma - 1))
   end
end

function flux_winters(equation::IsentropicEuler, lam, Ul, Ur, x, Uf) # Numerical flux of face at x
   gamma  = equation.γ
   @unpack epsilon = equation

    # Unpack left and right state
    rho_ll, v1_ll = pde2primitive(Ul, gamma)
    rho_rr, v1_rr = pde2primitive(Ur, gamma)
    kappa = 1.0

    p_ll = kappa * rho_ll^gamma
    p_rr = kappa * rho_rr^gamma

    # Compute the necessary mean values
    if gamma == 1.0 # isothermal gas
        rho_mean = ln_mean(rho_ll, rho_rr)
    else # equations.gamma > 1 # polytropic gas
        rho_mean = stolarsky_mean(rho_ll, rho_rr, gamma)
    end
    v1_avg = 0.5 * (v1_ll + v1_rr)
    p_avg = 0.5 * (p_ll + p_rr)

    f1 = rho_mean * 0.5 * (v1_ll + v1_rr)
    f2 = f1 * v1_avg + p_avg / epsilon^2

    return SVector(f1, f2)
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
                 "flux_winters"   => flux_winters,
                 "roe"            => roe!,
                 "hll"            => hll!,
                 "hllc"           => hllc!,
                 "flux_central"   => flux_central
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
   L = maximum(abs.(Up[1,1:nx]))
   p1 = @views plot(xc, Up[1,1:nx], legend=false, label = nothing,
                    ylim = (ymin-0.1*L, ymax+0.1*L),
                    linestyle = :dot, color = :red,
                    markerstrokestyle = :dot, markershape = :circle,
                    markersize = 3, markerstrokealpha = 0)
   xlabel!(p1, "x"); ylabel!(p1, "Density")

   ymin, ymax = minimum(Up[2,1:nx]), maximum(Up[2,1:nx])
   L = maximum(abs.(Up[2,1:nx]))
   p2 = @views plot(xc, Up[2,1:nx], legend=false, label = nothing,
                    ylim = (ymin-0.1*L, ymax+0.1*L),
                    linestyle = :dot, color = :red,
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
   # @assert false size(Up[:,1:nx])

   open("output/data_$time.txt", "w") do file
      write(file, "$t\n")
      writedlm(file, Up[:,1:nx]')
  end
   title!(p[1], "$nx points, time = $time")
   for i=1:nvar
      L = maximum(Up[i,1:nx]) - minimum(Up[i,1:nx])
      y_lims = (minimum(Up[i,:])-0.1*L, maximum(Up[i,:])+0.1*L)
      @show y_lims, L
      ylims!(p[i+1],y_lims) # Bad solution
      p[i+1][1][:y] = @views Up[i,1:nx]
   end
   frame(anim)
end

function final_plot(plt_data, equation::IsentropicEuler)
   p, anim = plt_data
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
get_equation(γ, epsilon) = (; eq = IsentropicEuler(γ, epsilon), flux, fprime, numfluxes,
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