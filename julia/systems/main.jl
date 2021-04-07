module main

push!(LOAD_PATH,".") # Temporarily done to ensure VTK_OUT works,
                     # will remove if not needed
using VTK_OUT
using LinearAlgebra
using OffsetArrays

# TODO - Add safety CFL
function compute_lam_dt(Ua, grid)
   nx, dx = grid.nx, grid.dx
   lam = 0.0
   dt  = 1.0
   for i=1:nx
      ua   = Ua[:, i]
      lam0 = abs.(eigvals(fprime(ua)))
      lam  = max(lam, lam0)
      dt   = min(dt, dx[i]/lam0)
   end
   return lam, dt
end

function set_initial_condition!(grid, U1)
   nx = grid.nx
   for i=1:nx
      U1[:,i] = initial_condition(xc[i])
   end
end

function set_initial_plot(grid,U1)
   xc = grid.xc
   @views plot(xc, U1[1,:], xc, U1[2,:], xc, U1[3,:])
end

gArray(nx, nvar) = OffsetArray(zeros(nx+2, nvar), 
                              OffsetArrays.Origin(0, 1))

function solve!(Grid, nvar)
   nx = Grid.nx
   dx = Grid.dx
   xf = Grid.xf
   # Allocating variables
   U1  = gArray(nvar, nx)
   res = gArray(nvar, nx) # dU/dt + res(U) = 0
   Ua  = U1 # ua is just Ua for this first order method,
            # storing for clarity
   set_initial_condition!(grid, U1)
   set_inial_plot(grid,U1)
   while t < Tf
      lam, dt = compute_lam_dt(Ua, grid)
      # loop over faces
      for i=1:nx+1
         Ul, Ur .= U1[:,i-1], U1[:,i]
         flux    = num_flux(Ul, Ur, eq, xf[i])
         @. res[:, i-1] += flux/dx[i]
         @. res[:, i]   -= flux/dx[i]
      end
      U1 -= dt*res
end