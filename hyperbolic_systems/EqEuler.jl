module EqEuler

using LinearAlgebra
using Plots
using DelimitedFiles
using LaTeXStrings

# For efficiency, euler should also contain γm1,γm3,3γm1_2. That just might a
# kill of readability though. Also, we don't know whether accessing that
# far-away gamma is even better than computing it every time
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

   A = [0.0                          1.0          0.0     ;
        0.5*(eq.γ-3.0)*u^2       (3.0-eq.γ)*u     eq.γ-1.0;
        u*(0.5*(eq.γ-1.0)*u^2-H) H-(eq.γ-1.0)*u^2 eq.γ*u  ]
   return A
end

# Numerical fluxes

function lax_friedrich(equation, lam, Ul, Ur, x) # Numerical flux of face at x
   eq = equation["eq"]
   Fl, Fr = flux(x, Ul, eq), flux(x, Ur, eq)
   value  = 0.5*(Fl+Fr) - 0.5 * lam * (Ur - Ul)
   return value
end

function steger_warming(equation, lam, Ul, Ur, x)
   eq = equation["eq"]
   γ  = eq.γ
   δ  = 0.0
   # Ul on Fp
   ρ, u, E = Ul[1], Ul[2]/Ul[1], Ul[3]    # density, velocity, energy
   p = (γ - 1.0) * (E - 0.5*ρ*u^2)        # pressure
   c = sqrt(γ*p/ρ)                        # sound speed
   λ = [u, u+c, u-c] # Poor ordering, should fix.
   λp = 0.5 * (λ + sqrt.(λ.^2 .+ δ^2))    # positive part of eigenvalues
   w  = 0.5*(3.0-γ)*(λp[2]+λp[3])*c^2/(γ-1.0)
   Fp = 0.5*ρ/γ * [ 2.0*(γ-1.0)*λp[1]   +     λp[2]         +     λp[3]
                    2.0*(γ-1.0)*λp[1]*u +     λp[2]*(u+c)   +     λp[3]*(u-c)
                    (γ-1.0)*λp[1]*u^2   + 0.5*λp[2]*(u+c)^2 + 0.5*λp[3]*(u-c)^2 + w  ]
   # Ur on Fm
   ρ, u, E = Ur[1], Ur[2]/Ur[1], Ur[3]    # density, velocity, energy
   p = (γ - 1.0) * (E - 0.5*ρ*u^2)        # pressure
   c = sqrt(γ*p/ρ)                        # sound speed
   λ = [u, u+c, u-c]                      # Poor ordering, should fix.
   λm = 0.5 * (λ - sqrt.(λ.^2 .+ δ^2))    # negative part of eigenvalues
   w = 0.5*(3.0-γ)*(λm[2]+λm[3])*c^2/(γ-1.0)
   Fm = 0.5*ρ/γ * [ 2.0*(γ-1.0)*λm[1]       +     λm[2]         +     λm[3]
                    2.0*(γ-1.0)*λm[1]*u     +     λm[2]*(u+c)   +     λm[3]*(u-c)
                        (γ-1.0)*λm[1]*u^2   + 0.5*λm[2]*(u+c)^2 + 0.5*λm[3]*(u-c)^2 + w ]
   # Should we make a subfunction to avoid duplication, or would it cause
   # performance issues?
   return Fp + Fm
end

# function converting primitive variables to PDE variables
function primitive2pde(prim, γ) # primitive, viscosity
   U = [prim[1], prim[1]*prim[2], prim[3]/(γ-1.0) + prim[1]*prim[2]^2/2.0]
      # ρ    ,     ρ*u     ,        p/(γ-1.0) +     ρ*u^2/2.0
   return U
end

# function converting pde variables to primitive variables
function pde2primitive(U, γ)
   primitives = [U[1], U[2]/U[1], (γ-1.0)*(U[3]-U[2]^2/(2.0*U[1]))]
   #            [ρ ,   u        , p]
   return primitives
end

function initialize_plot(grid, problem, equation, U)
   anim = Animation()
   xc = grid.xc
   nx = grid.nx
   eq = equation["eq"]
   nvar = problem["nvar"]
   Up = copy(U)
   for j=1:nx
      @views Up[:, j] = pde2primitive(U[:,j],eq.γ)
   end
   # Adding title as a subplot in itself
   p_title = title = plot(title = "Solution at time = 0", grid = false,
                          showaxis = false, bottom_margin = -50Plots.px)
   ymin, ymax = minimum(Up[1,1:nx]), maximum(Up[1,1:nx])
   p1 = @views plot(xc, Up[1,1:nx], legend=false, label = nothing,
                    ylim = (ymin-0.1, ymax+0.1))
   title!(p1, "Density")
   xlabel!(p1, L"x"); ylabel!(p1, L"U")

   ymin, ymax = minimum(Up[2,1:nx]), maximum(Up[2,1:nx])
   p2 = @views plot(xc, Up[2,1:nx], legend=false, label = nothing,
                    ylim = (ymin-0.1, ymax+0.1))
   title!(p2, "Velocity")
   xlabel!(p2, L"x"); ylabel!(p2, L"U")

   ymin, ymax = minimum(Up[3,1:nx]), maximum(Up[3,1:nx])
   p3 = @views plot(xc, Up[3,1:nx], legend=false, label = nothing,
                    ylim = (ymin-0.1, ymax+0.1))
   title!(p3, "Pressure")
   xlabel!(p3, L"x"); ylabel!(p3, L"U")

   l = @layout[ a{0.01h}; b c d]
   p = plot(p_title, p1, p2, p3, layout = l) # Can this be nvar independent?
   frame(anim)
   return p, anim
end

function update_plot!(grid, equation, problem, U, t, it, param, p, anim)
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
   eq = equation["eq"]
   nvar = problem["nvar"]
   Up = copy(U)
   for j=1:nx
      @views Up[:, j] = pde2primitive(U[:,j],eq.γ)
   end
   title!(p[1], "Solution at time "*string(t))
   for i=1:nvar
      y_lims = (minimum(Up[i,:])-0.1, maximum(Up[i,:])+0.1)
      ylims!(p[i+1],y_lims) # Bad solution
      p[i+1][1][:y] = @views Up[i,1:nx]
   end
   frame(anim)
end

get_equation(γ) = Dict( "eq"              => Euler(γ),
                        "flux"            => flux,
                        "fprime"          => fprime,
                        "initialize_plot" => initialize_plot,
                        "update_plot!"    => update_plot!,
                        "name"            => "1D Euler equations")


export lax_friedrich
export steger_warming
export get_equation
export primitive2pde
export pde2primitive
export plot_final_soln

end
