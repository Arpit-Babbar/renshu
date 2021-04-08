module main

using LinearAlgebra
using OffsetArrays
using PyPlot

gArray(nvar, nx) = OffsetArray(zeros(nvar, nx+2), 
                              OffsetArrays.Origin(1, 0))

# TODO - Add safety CFL
function compute_lam_dt(eq, grid, Ua)
   nx, dx = grid.nx, grid.dx
   xc     = grid.xc
   lam = 0.0
   dt  = 1.0
   for i=1:nx
      ua   = Ua[:, i]
      lam0 = maximum(abs.(eigvals(eq.fprime(ua, xc[i])))) # CHECK xc[i]!
      lam  = max(lam, lam0)
      dt   = min(dt, dx[i]/lam0)
   end
   println("dt = ", dt)
   return lam, 0.5*dt
end

function set_initial_condition!(grid, U, initial_condition)
   nx = grid.nx
   xc = grid.xc
   for i=1:nx
      U[:,i] = initial_condition(xc[i])
   end
end

function plot_solution(grid, U, Ue, t)
   plt. clf()
   xc = grid.xc
   nx = grid.nx
   title("Solutions plot at time $t")
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

function compute_exact_soln!(eq, grid, t, initial_condition, nvar, Ue)
   nx = grid.nx
   xc = grid.xc

   fp = eq.fprime(Ue, 1.0) # 1.0 is dummy
   eigen_decomp = eigen(fp)
   lam, eigen_vecs = eigen_decomp.values, eigen_decomp.vectors
   for j=1:nx
      for i=1:nvar
         # Too complicated, fix.
         Ue[i,j] = (inv(eigen_vecs) * initial_condition(xc[j] - lam[i] * t))[i]
      end
   end
   for j=1:nx
      Ue[:,j] = eigen_vecs * Ue[:,j]
   end
   return nothing
end

function update_ghost!(grid, U, Ue)
   nx = grid.nx
   U[:, 0]    .= Ue[:, 1]
   U[:, nx+1] .= Ue[:, nx]
   return nothing
end

function compute_residual!(eq, grid, lam, U, num_flux, res)
   nx = grid.nx
   xf = grid.xf
   dx = grid.dx
   dx0 =  OffsetArray(zeros(nx+2), OffsetArrays.Origin(0))
   dx0[1:nx] .= dx
   # loop over faces
   for i=1:nx+1
      Ul, Ur  = U[:,i-1], U[:,i]
      f    = num_flux(eq, lam, Ul, Ur, xf[i])
      @views res[:, i-1] += f/dx0[i-1]
      @views res[:, i]   -= f/dx0[i]
   end
end

function solve(eq, grid, initial_condition, num_flux, Tf, nvar)
   nx = grid.nx
   dx = grid.dx
   xf = grid.xf
   # Allocating variables
   U   = gArray(nvar, nx)
   Ue  = zeros(nvar, nx)
   res = gArray(nvar, nx) # dU/dt + res(U) = 0
   Ua  = U # ua is just Ua for this first order method,
            # storing for clarity
   set_initial_condition!(grid, U, initial_condition)
   it, t = 0, 0.0
   figure(figsize=(15,5))
   while t < Tf
      # lam, dt = compute_lam_dt(eq, grid, Ua)
      lam, dt = 3.0, dx[1]/10.0
      compute_exact_soln!(eq, grid, t, initial_condition, nvar, Ue)
      update_ghost!(grid, U, Ue)                            # Fills ghost cells
      compute_residual!(eq, grid, lam, U, num_flux, res)
      @. U -= dt*res
      plot_solution(grid, U, Ue, t)
      t += dt; it += 1
   end
end

export solve

end