module EqEuler

using LinearAlgebra
using Plots
using DelimitedFiles
using StaticArrays
using LaTeXStrings
using UnPack

import ..FV1D.FV: get_node_vars

# For efficiency, euler should also contain γm1,γm3,3γm1_2. That just might a
# kill of readability though. Also, we don't know whether accessing that
# far-away gamma is even better than computing it every time
struct Euler
   γ::Float64
end

# Would have liked to solver

id_(x) = x

# pure function
function flux(x, U, eq::Euler)
   γ = eq.γ
   ρ = U[1]        # density
   u = U[2] / U[1] # velocity
   E = U[3]        # energy
   p = (γ - 1.0) * (E - 0.5 * ρ * u^2) # pressure
   F = SVector(U[2], p + ρ * u^2, (E+p) * u) # flux
   return F
end

function get_node_vars(U, eq::Euler, indices)
   SVector(ntuple(v -> U[v, indices], 3))
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

   U = SVector(prim[1], prim[1]*prim[2], prim[3]/(γ-1.0) + prim[1]*prim[2]^2/2.0)
      # ρ    ,     ρ*u     ,        p/(γ-1.0) +     ρ*u^2/2.0
   return U
end

# function converting pde variables to primitive variables
function pde2primitive(U, γ)
   primitives = SVector(U[1], U[2]/U[1], (γ-1.0)*(U[3]-U[2]^2/(2.0*U[1])))
   #            [ρ ,   u        , p]
   return primitives
end

#-------------------------------------------------------------------------------
# Numerical Fluxes
#-------------------------------------------------------------------------------
function lax_friedrich!(equation::Euler, lam, Ul, Ur, x, Uf) # Numerical flux of face at x
   @unpack eq = equation
   Fl, Fr = flux(x, Ul, eq), flux(x, Ur, eq)
   Uf  .= 0.5*(Fl+Fr) - 0.5 * lam * (Ur - Ul)
   return nothing
end

function rusanov!(equation::Euler, lam, Ul, Ur, x, Uf) # Numerical flux of face at x
   γ  = equation.γ
   ρl, ul, pl = pde2primitive(Ul, γ)
   ρr, ur, pr = pde2primitive(Ur, γ)
   cl, cr = sqrt(γ*pl/ρl), sqrt(γ*pr/ρr)                   # sound speed
   λ = maximum(abs.([ul, ul-cl, ul+cl, ur, ur-cr, ur+cr])) # local wave speed
   Fl, Fr = flux(x, Ul, equation), flux(x, Ur, equation)
   F  = 0.5*(Fl+Fr) - 0.5*λ*(Ur - Ul)
   return F
end

function steger_warming!(equation::Euler, lam, Ul, Ur, x, Uf)
   @unpack eq = equation
   γ  = eq.γ
   δ  = 0.0
   # Ul on Fp
   ρ, u, E = Ul[1], Ul[2]/Ul[1], Ul[3]    # density, velocity, energy
   p = (γ - 1.0) * (E - 0.5*ρ*u^2)        # pressure
   c = sqrt(γ*p/ρ)                        # sound speed
   λ = [u, u+c, u-c] # Poor ordering, should fix.
   λp = 0.5 * (λ + sqrt.(λ.^2 .+ δ^2))    # positive part of eigenvalues
   w  = 0.5*(3.0-γ)*(λp[2]+λp[3])*c^2/(γ-1.0)
   Uf .= 0.5*ρ/γ * [ 2.0*(γ-1.0)*λp[1]   +     λp[2]         +     λp[3]
                    2.0*(γ-1.0)*λp[1]*u +     λp[2]*(u+c)   +     λp[3]*(u-c)
                    (γ-1.0)*λp[1]*u^2   + 0.5*λp[2]*(u+c)^2 + 0.5*λp[3]*(u-c)^2 + w ]
   # Uf  = Fp
   # Ur on Fm
   ρ, u, E = Ur[1], Ur[2]/Ur[1], Ur[3]    # density, velocity, energy
   p = (γ - 1.0) * (E - 0.5*ρ*u^2)        # pressure
   c = sqrt(γ*p/ρ)                        # sound speed
   λ = [u, u+c, u-c]                      # Poor ordering, should fix.
   λm = 0.5 * (λ - sqrt.(λ.^2 .+ δ^2))    # negative part of eigenvalues
   w = 0.5*(3.0-γ)*(λm[2]+λm[3])*c^2/(γ-1.0)
   Uf .+= 0.5*ρ/γ * [ 2.0*(γ-1.0)*λm[1]       +     λm[2]         +     λm[3]
                    2.0*(γ-1.0)*λm[1]*u     +     λm[2]*(u+c)   +     λm[3]*(u-c)
                        (γ-1.0)*λm[1]*u^2   + 0.5*λm[2]*(u+c)^2 + 0.5*λm[3]*(u-c)^2 + w ]
   # Uf = Fp + Fm
   # Should we make a subfunction to avoid duplication, or would it cause
   # performance issues?
   return nothing
end

function roe!(equation::Euler, lam, Ul, Ur, x, Uf)
   @unpack eq = equation
   γ  = eq.γ
   ϵ  = 0.2
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
   Uf .= 0.5*(Fl+Fr) - 0.5*(α1*l1*r1 + α2*l2*r2 + α3*l3*r3)
   return nothing
end

function vanleer!(equation::Euler, lam, Ul, Ur, x, Uf)
   @unpack eq = equation
   γ  = eq.γ
   # Ul on Fp
   ρ, u, E = Ul[1], Ul[2]/Ul[1], Ul[3]    # density, velocity, energy
   p = (γ - 1.0) * (E - 0.5*ρ*u^2)        # pressure
   a = sqrt(γ*p/ρ)                        # sound speed
   M = u/a                                # mach number
   Fl = flux(x, Ul, eq)
   if M <= -1
      # Uf .= zeros(Float64, 3)
      nothing
   elseif M >= 1
      Uf .= Fl
      return nothing
   else
      Uf .= 0.25*ρ*a*(1.0+M)^2 * [1.0
                                 2.0*a/γ           * (0.5*(γ-1.0)*M+1.0)
                                 2.0*a^2/(γ^2-1.0) * (0.5*(γ-1.0)*M+1.0)^2]
      # Uf = Fp
   end
   # Ur on Fm
   ρ, u, E = Ur[1], Ur[2]/Ur[1], Ur[3]    # density, velocity, energy
   p = (γ - 1.0) * (E - 0.5*ρ*u^2)        # pressure
   a = sqrt(γ*p/ρ)                        # sound speed
   M = u/a                                # mach number
   Fr = flux(x, Ur, eq)
   if M <= -1
      Uf .+= Fr
   elseif M >= 1
      nothing
      # Uf .+= zeros(Float64, 3)
   else
      Uf .+= -0.25*ρ*a*(1.0-M)^2 * [1.0
                                 2.0*a/γ           * (0.5*(γ-1.0)*M-1.0)
                                 2.0*a^2/(γ^2-1.0) * (0.5*(γ-1.0)*M-1.0)^2]
   end
   return nothing
end

function hll!(equation::Euler, lam, Ul, Ur, x, Uf)
   γ = equation.γ
   # TODO - Replace with pde2primitive
   ρl, ul, El = Ul[1], Ul[2]/Ul[1], Ul[3]    # density, velocity, energy
   pl = (γ - 1.0) * (El - 0.5*ρl*ul^2)        # pressure
   cl = sqrt(γ*pl/ρl)                        # sound speed
   ρr, ur, Er = Ur[1], Ur[2]/Ur[1], Ur[3]    # density, velocity, energy
   pr = (γ - 1.0) * (Er - 0.5*ρr*ur^2)        # pressure
   cr = sqrt(γ*pr/ρr)                        # sound speed
   Hl, Hr = γ*pl/((γ-1.0)*ρl) + 0.5*ul^2 , γ*pr/((γ-1.0)*ρr) + 0.5*ur^2 # enthl
   ⎷ρl, ⎷ρr = sqrt(ρl), sqrt(ρr) # for efficiency
   u = (⎷ρl*ul + ⎷ρr*ur) / (⎷ρl + ⎷ρr)      # roe avg velocity
   H = (⎷ρl*Hl + ⎷ρr*Hr) / (⎷ρl + ⎷ρr)      # roe avg enthalpy
   c = sqrt((γ-1.0) * (H - 0.5*u^2))         # sound speed
   Sl, Sr = min(ul-cl,u-c), max(ur+cr,u+c)
   Fl, Fr = flux(x, Ul, equation), flux(x, Ur, equation)
   dU = Ur - Ul
   dS = Sr - Sl
   if Sl > 0
      Uf .= Fl
   elseif Sr < 0
      Uf .= Fr
   else
      Uf .= (Sr*Fl-Sl*Fr + Sl*Sr*dU)/dS
   end
   return nothing
end


function hllc!(equation::Euler, lam, Ul, Ur, x, Uf)
   @unpack eq = equation
   γ = eq.γ
   # TODO - Replace with pde2primitive
   # Choice of Sl, Sr from Einfeldt et al.
   ρl, ul, El = Ul[1], Ul[2]/Ul[1], Ul[3]    # density, velocity, energy
   pl = (γ - 1.0) * (El - 0.5*ρl*ul^2)        # pressure
   cl = sqrt(γ*pl/ρl)                        # sound speed
   ρr, ur, Er = Ur[1], Ur[2]/Ur[1], Ur[3]    # density, velocity, energy
   pr = (γ - 1.0) * (Er - 0.5*ρr*ur^2)        # pressure
   cr = sqrt(γ*pr/ρr)                        # sound speed
   Hl, Hr = γ*pl/((γ-1.0)*ρl) + 0.5*ul^2 , γ*pr/((γ-1.0)*ρr) + 0.5*ur^2 # enthl
   ⎷ρl, ⎷ρr = sqrt(ρl), sqrt(ρr) # for efficiency
   u = (⎷ρl*ul + ⎷ρr*ur) / (⎷ρl + ⎷ρr)      # roe avg velocity
   H = (⎷ρl*Hl + ⎷ρr*Hr) / (⎷ρl + ⎷ρr)      # roe avg enthalpy
   c = sqrt((γ-1.0) * (H - 0.5*u^2))         # sound speed
   Sl, Sr = min(ul-cl,u-c), max(ur+cr,u+c)
   Fl, Fr = flux(x, Ul, eq), flux(x, Ur, eq)
   if Sl >= 0.0
      Uf .= Fl
      return nothing
   elseif Sr <= 0.0
      Uf .= Fr
      return nothing
   end
   # u✶
   Δp = pr-pl
   ustar = (ρr*ur*(Sr-ur)-ρl*ul*(Sl-ul)-Δp)/(ρr*(Sr-ur)-ρl*(Sl-ul))
   Sstar = ustar
   # ρstar_l, ρstar_r
   ρstar_l = (Sl-ul)/(Sl-Sstar) * ρl
   ρstar_r = (Sr-ur)/(Sr-Sstar) * ρr
   # pstar
   pstar = pl + ρl*(Sl-ul)*(ustar-ul)
   # pstar = 0.5*(pl+pr) + 0.5*(ρl*(Sl-ul)*(ustar-ul)+ρr*(Sr-ur)*(ustar-ur))
   # Estar_l, Estar_r
   # Estar_l = ρstar_l*(El/ρl+(Sstar-ul)*(Sstar+pl/(ρl*(Sl-ul))))
   # Estar_r = ρstar_r*(Er/ρr+(Sstar-ur)*(Sstar+pr/(ρr*(Sr-ur))))
   Estar_l = ((Sl-ul)*El+pstar*ustar-pl*ul)/(Sl-ustar)
   Estar_r = ((Sr-ur)*Er+pstar*ustar-pr*ur)/(Sr-ustar)
   Ustar_l = [ρstar_l, ρstar_l*ustar, Estar_l]
   Ustar_r = [ρstar_r, ρstar_r*ustar, Estar_r]
   if ustar >= 0.0
      Uf .= Fl+Sl*(Ustar_l-Ul)
      return nothing
   elseif ustar < 0.0
      Uf .= Fr+Sr*(Ustar_r-Ur)
      return nothing
   end
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
function initialize_plot(grid, problem, equation::Euler, scheme, U)
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

function update_plot!(grid, problem, equation::Euler, scheme, U, t, it, param, plt_data)
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

function final_plot(plt_data, equation::Euler)
   p, anim = plt_data
   soln_data = readdlm("toro_user_exact.dat", skipstart = 9);
   @views x = soln_data[:,1];
   @views dens_exact = soln_data[:,2];
   @views pres_exact = soln_data[:,3];
   @views velx_exact = soln_data[:,4];
   plot!(p[2],x,dens_exact, label = nothing, color = :blue, legend=false)
   plot!(p[4],x,pres_exact, label = nothing, color = :blue, legend=false)
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

# get_equation(γ) = Dict( "eq"              => Euler(γ),
#                         "flux"            => flux,
#                         "fprime"          => fprime,
#                         "numfluxes"       => numfluxes,
#                         "name"            => "1D Euler equations")
get_equation(γ) = (; eq = Euler(γ), flux, fprime, numfluxes,
                     name = "1D Euler equations")

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