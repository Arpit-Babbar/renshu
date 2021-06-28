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

# Numerical fluxes

function lax_friedrich(equation, lam, Ul, Ur, x) # Numerical flux of face at x
   eq = equation["eq"]
   Fl, Fr = flux(x, Ul, eq), flux(x, Ur, eq)
   value  = 0.5*(Fl+Fr) - 0.5 * lam * (Ur - Ul)
   return value
end

function rusanov(equation, lam, Ul, Ur, x) # Numerical flux of face at x
   eq = equation["eq"]
   γ  = eq.γ
   ρl, ul, pl = pde2primitive(Ul, γ)
   ρr, ur, pr = pde2primitive(Ur, γ)
   cl, cr = sqrt(γ*pl/ρl), sqrt(γ*pr/ρr)                   # sound speed
   λ = maximum(abs.([ul, ul-cl, ul+cl, ur, ur-cr, ur+cr])) # local wave speed
   Fl, Fr = flux(x, Ul, eq), flux(x, Ur, eq)
   value  = 0.5*(Fl+Fr) - 0.5*λ*(Ur - Ul)
   return value
end

# Add Rusanov as well!

function steger_warming(equation, lam, Ul, Ur, x)
   eq = equation["eq"]
   γ  = eq.γ
   δ  = 0.1
   # Ul on Fp
   ρ, u, E = Ul[1], Ul[2]/Ul[1], Ul[3]    # density, velocity, energy
   p = (γ - 1.0) * (E - 0.5*ρ*u^2)        # pressure
   c = sqrt(γ*p/ρ)                        # sound speed
   λ = [u, u+c, u-c] # Poor ordering, should fix.
   λp = 0.5 * (λ + sqrt.(λ.^2 .+ δ^2))    # positive part of eigenvalues
   w  = 0.5*(3.0-γ)*(λp[2]+λp[3])*c^2/(γ-1.0)
   Fp = 0.5*ρ/γ * [ 2.0*(γ-1.0)*λp[1]   +     λp[2]         +     λp[3]
                    2.0*(γ-1.0)*λp[1]*u +     λp[2]*(u+c)   +     λp[3]*(u-c)
                    (γ-1.0)*λp[1]*u^2   + 0.5*λp[2]*(u+c)^2 + 0.5*λp[3]*(u-c)^2 + w ]
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
   output = Fp + Fm
   return output
end

function roe(equation, lam, Ul, Ur, x)
   eq = equation["eq"]
   γ  = eq.γ
   ϵ  = 0.0
   ρl, ul, El = Ul[1], Ul[2]/Ul[1], Ul[3]    # density, velocity, energy
   ρr, ur, Er = Ur[1], Ur[2]/Ur[1], Ur[3]    # density, velocity, energy
   pl, pr  = (γ - 1.0)*(El - 0.5*ρl*ul^2), (γ - 1.0)*(Er - 0.5*ρr*ur^2) # press
   ⎷ρl, ⎷ρr = sqrt(ρl), sqrt(ρr) # for efficiency
   Hl, Hr = γ*pl/((γ-1.0)*ρl) + 0.5*ul^2 , γ*pr/((γ-1.0)*ρr) + 0.5*ur^2 # enthl
   u = (⎷ρl*ul + ⎷ρr*ur) / (⎷ρl + ⎷ρr)      # roe avg velocity
   H = (⎷ρl*Hl + ⎷ρr*Hr) / (⎷ρl + ⎷ρr)      # roe avg enthalpy
   c = sqrt((γ-1.0) * (H - 0.5*u^2))         # sound speed
   # Eigenvectors
   r1 = [1.0, u - c, H - u*c ]
   r2 = [1.0, u,     0.5*u^2 ]
   r3 = [1.0, u + c, H + u*c ]
   # Computing R |L| inv(R) ΔU efficiently
   dU = Ur - Ul
   α2 = (γ - 1.0)/c^2 * ( (H - u^2)*dU[1] + u*dU[2] - dU[3] )
   α1 = 1.0/(2.0 * c) * ( (u+c) * dU[1] - dU[2] - c*α2 )
   α3 = dU[1] - α1 - α2
   l1, l2, l3 = abs(u-c), abs(u), abs(u+c)
   δ  = c*ϵ
   if abs(l1)<2.0*ϵ l1 = l1^2/(4.0*δ) + δ end
   if abs(l2)<2.0*ϵ l2 = l2^2/(4.0*δ) + δ end
   if abs(l3)<2.0*ϵ l3 = l3^2/(4.0*δ) + δ end
   # compute flux
   Fl, Fr = flux(x, Ul, eq), flux(x, Ur, eq)
   output = 0.5*(Fl+Fr) - 0.5*(α1*l1*r1 + α2*l2*r2 + α3*l3*r3)
   return output
end

function vanleer(equation, lam, Ul, Ur, x)
   eq = equation["eq"]
   γ  = eq.γ
   # Ul on Fp
   ρ, u, E = Ul[1], Ul[2]/Ul[1], Ul[3]    # density, velocity, energy
   p = (γ - 1.0) * (E - 0.5*ρ*u^2)        # pressure
   a = sqrt(γ*p/ρ)                        # sound speed
   M = u/a                                # mach number
   Fp = 0.25*ρ*a*(1.0+M)^2 * [1.0
                              2.0*a/γ           * (0.5*(γ-1.0)*M+1.0)
                              2.0*a^2/(γ^2-1.0) * (0.5*(γ-1.0)*M+1.0)^2]
   # Ur on Fm
   ρ, u, E = Ur[1], Ur[2]/Ur[1], Ur[3]    # density, velocity, energy
   p = (γ - 1.0) * (E - 0.5*ρ*u^2)        # pressure
   a = sqrt(γ*p/ρ)                        # sound speed
   M = u/a                                # mach number
   Fm = -0.25*ρ*a*(1.0-M)^2 * [1.0
                               2.0*a/γ           * (0.5*(γ-1.0)*M-1.0)
                               2.0*a^2/(γ^2-1.0) * (0.5*(γ-1.0)*M-1.0)^2]
   output = Fp + Fm
   return output
end

function hl(equation, lam, Ul, Ur, x)
   eq = equation["eq"]
   γ = eq.γ
   # TODO - Replace with pde2primitive
   ρl, ul, El = Ul[1], Ul[2]/Ul[1], Ul[3]    # density, velocity, energy
   pl = (γ - 1.0) * (El - 0.5*ρl*ul^2)        # pressure
   cl = sqrt(γ*pl/ρl)                        # sound speed
   ρr, ur, Er = Ul[1], Ul[2]/Ul[1], Ul[3]    # density, velocity, energy
   pr = (γ - 1.0) * (Er - 0.5*ρr*ur^2)        # pressure
   cr = sqrt(γ*pr/ρr)                        # sound speed
   ⎷ρl, ⎷ρr = sqrt(ρl), sqrt(ρr) # for efficiency
   u = (⎷ρl*ul + ⎷ρr*ur) / (⎷ρl + ⎷ρr)      # roe avg velocity
   p = (⎷ρl*pl + ⎷ρr*pr) / (⎷ρl + ⎷ρr)      # roe avg pressure
   ρ = (⎷ρl*ρl + ⎷ρr*ρr) / (⎷ρl + ⎷ρr)      # roe avg pressure
   c = sqrt(γ*p/ρ)                           # roe avg speed
   Sl, Sr = min(ul-cl,u-c), max(ur+cr,u+c)
   Fl, Fr = flux(x, Ul, eq), flux(x, Ur, eq)
   dU = Ur - Ul
   dS = Sr - Sl
   output = (Sr*Fl-Sl*Fr + Sl*Sr*dU)/dS
   return output
end

function initialize_plot(grid, problem, equation, scheme, U)
   anim = Animation()
   xc = grid.xc
   nx = grid.nx
   eq = equation["eq"]
   nvar = problem["nvar"]
   numflux = scheme["numflux_ind"]      # numflux as string
   numflux = replace(numflux, "_"=>" ") # Remove underscore
   numflux = titlecase(numflux)         # Capitalize
   Up = copy(U)
   for j=1:nx
      @views Up[:, j] = pde2primitive(U[:,j],eq.γ)
   end
   # Adding title as a subplot in itself
   p_title = plot(title = "$numflux flux, $nx points, time = 0",
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

   ymin, ymax = minimum(Up[3,1:nx]), maximum(Up[3,1:nx])
   p3 = @views plot(xc, Up[3,1:nx], legend=true, label = "Approximate",
                    ylim = (ymin-0.1, ymax+0.1), linestyle = :dot, color = :red,
                    markerstrokestyle = :dot, markershape = :circle,
                    markersize = 3, markerstrokealpha = 0)
   xlabel!(p3, "x"); ylabel!(p3, "Pressure")

   l = @layout[ a{0.01h}; b c d]
   p = plot(p_title, p1, p2, p3, layout = l) # Can this be nvar independent?
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
   eq = equation["eq"]
   nvar = problem["nvar"]
   numflux = scheme["numflux_ind"]      # numflux as string
   numflux = replace(numflux, "_"=>" ") # Remove underscore
   numflux = titlecase(numflux)         # Capitalize
   Up = copy(U)
   for j=1:nx
      @views Up[:, j] = pde2primitive(U[:,j],eq.γ)
   end
   time = round(t, digits=3)
   title!(p[1], "$numflux flux, $nx points, time = $time")
   for i=1:nvar
      y_lims = (minimum(Up[i,:])-0.1, maximum(Up[i,:])+0.1)
      ylims!(p[i+1],y_lims) # Bad solution
      p[i+1][1][:y] = @views Up[i,1:nx]
   end
   frame(anim)
end

numfluxes = Dict("lax_friedrich"  => lax_friedrich,
                 "rusanov"        => rusanov,
                 "steger_warming" => steger_warming,
                 "vanleer"        => vanleer,
                 "roe"            => roe,
                 "hl"             => hl
                 )

get_equation(γ) = Dict( "eq"              => Euler(γ),
                        "flux"            => flux,
                        "fprime"          => fprime,
                        "initialize_plot" => initialize_plot,
                        "update_plot!"    => update_plot!,
                        "numfluxes"       => numfluxes,
                        "name"            => "1D Euler equations")

export roe
export lax_friedrich
export hl
export rusanov
export steger_warming
export vanleer
export get_equation
export primitive2pde
export pde2primitive
export plot_final_soln

end
