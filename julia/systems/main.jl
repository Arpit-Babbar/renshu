module main

push!(LOAD_PATH,".") # Temporarily done to ensure VTK_OUT works,
                     # will remove if not needed
using VTK_OUT
using LinearAlgebra
using OffsetArrays

gArray(nx, nvar) = OffsetArray(zeros(nx+2, nvar), 
                              OffsetArrays.Origin(0, 1))

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

function set_initial_condition!(grid, U, initial_condition)
   nx = grid.nx
   for i=1:nx
      U[:,i] = initial_condition(xc[i])
   end
end

function plot_solution(grid, U, t)
   xc = grid.xc
   figure(figsize=(15,5))
   title("Solutions plot at time $t")
   subplot(131)
   @views plot(xc, U[1,:])
   xlabel("x")
   ylabel("\$U_1\$")
   subplot(132)
   @views plot(xc, U[2,:])
   xlabel("x")
   ylabel("\$U_2\$")
   subplot(133)
   @views plot(xc, U[3,:])
   xlabel("x")
   ylabel("\$U_3\$")
end

function compute_exact_soln!(eq, grid, t, initial_condition, Ue)
   xc = grid.xc
   fp = eq.fprime(Ue, 1.0) # 1.0 is dummy
   eigen_decomp = eigen(fp)
   lam, eigen_vecs = eigen_decomp.values, eigen_decomp.vectors
   Ue = inv(eigen_vecs) * initial_condition(xc .- lam[i] * t)
   return nothing
end

function update_ghost!(grid, U, Ue)
   nx = grid.nx
   U[0]    = Ue[1]
   U[nx+1] = Ue[nx]
   return nothing
end

function compute_residual!(grid, U, res)
   nx = grid.nx
   # loop over faces
   for i=1:nx+1
      Ul, Ur .= U[:,i-1], U[:,i]
      flux    = num_flux(Ul, Ur, eq, xf[i])
      @views res[:, i-1] += flux/dx[i]
      @views res[:, i]   -= flux/dx[i]
   end
end

function solve!(eq, Grid, nvar)
   nx = Grid.nx
   dx = Grid.dx
   xf = Grid.xf
   # Allocating variables
   U   = gArray(nvar, nx)
   Ue  = zeros(nvar, nx)
   res = gArray(nvar, nx) # dU/dt + res(U) = 0
   Ua  = U # ua is just Ua for this first order method,
            # storing for clarity
   set_initial_condition!(grid, U, initial_condition)
                                                          # using exact solution
   set_inial_plot(grid,U)
   it, t = 0, 0.0
   while t < Tf
      lam, dt = compute_lam_dt(Ua, grid)
      compute_exact_soln!(eq, grid, t, initial_condition, Ue)
      update_ghost!(grid, U, Ue)                             # Fills ghost cells
      plot_solution(grid, U,t)
      @. U -= dt*res
      t += dt; iter += 1
end