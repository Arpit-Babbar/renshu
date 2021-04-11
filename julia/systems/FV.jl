module FV

using Grid
using LinearAlgebra
using OffsetArrays
using PyPlot

#-------------------------------------------------------------------------------
# Create a dictionary of problem description
#-------------------------------------------------------------------------------
function Problem(domain::Tuple{Float64,Float64},
                 nvar::Int64,
                 initial_value::Function,
                 boundary_value::Function,
                 boundary_condition::Any,
                 final_time::Float64)
Dict("domain" => domain,
     "nvar"   => nvar,
     "initial_value" => initial_value,
     "boundary_value" => boundary_value,
     "boundary_condition" => boundary_condition,
     "final_time" => final_time)
end

#-------------------------------------------------------------------------------
# Create a dictionary of parameters
#-------------------------------------------------------------------------------
function Parameters(grid_size::Int64,
                    cfl::Float64,
                    save_time_interval::Float64)
@assert (cfl >= 0.0) "cfl must be >= 0.0"
@assert (save_time_interval >= 0.0) "save_time_interval must be >= 0.0"
Dict("grid_size" => grid_size,
     "cfl" => cfl,
     "save_time_interval" => save_time_interval)
end

#-------------------------------------------------------------------------------
# Create a dictionary of scheme description
#-------------------------------------------------------------------------------
function Scheme(numerical_flux::Function)
   Dict("numerical_flux" => numerical_flux)
end

#-------------------------------------------------------------------------------
# Create ghosted arrays with one layer of cells on both sides
#-------------------------------------------------------------------------------
gArray(nvar, nx) = OffsetArray(zeros(nvar, nx+2), 
                              OffsetArrays.Origin(1, 0))

#-------------------------------------------------------------------------------
# Adjust dt to reach final time or the next time when solution has to be saved
#-------------------------------------------------------------------------------
function adjust_time_step(problem, param, dt, t)
   # Adjust to reach final time exactly
   final_time = problem["final_time"]
   save_time_interval = param["save_time_interval"]
   if t + dt > final_time
      dt = final_time - t
      return dt
   end

   # Adjust to reach next solution saving time
   if save_time_interval > 0.0
      next_save_time = ceil(t/save_time_interval) * save_time_interval
      # If t is not a plotting time, we check if the next time
      # would step over the plotting time to adjust dt
      if abs(t-next_save_time) > 1e-10 && t + dt - next_save_time > -1e-10
         dt = next_save_time - t
         return dt
      end
   end

   return dt
end

# TODO - Add safety CFL
function compute_lam_dt(equation, grid, Ua)
   eq = equation["eq"]
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
   return lam, dt
end

function set_initial_value!(grid, U, problem)
   nx = grid.nx
   xc = grid.xc
   initial_value = problem["initial_value"]
   for i=1:nx
      U[:,i] = initial_value(xc[i])
   end
end

function plot_solution(grid, U, Ue, t, it, param)
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

function update_ghost!(grid, U, Ue)
   nx = grid.nx
   U[:, 0]    .= Ue[:, 1]
   U[:, nx+1] .= Ue[:, nx]
   return nothing
end

function compute_residual!(equation, grid, lam, U, scheme, res)
   nx = grid.nx
   xf = grid.xf
   dx = grid.dx
   eq = equation["eq"]
   num_flux = scheme["numerical_flux"]
   dx0 =  OffsetArray(zeros(nx+2), OffsetArrays.Origin(0))
   dx0[1:nx] .= dx
   dx0[0] = dx[nx]
   dx0[nx+1] = dx[1]
   res[:,:] .= 0.0
   # loop over faces
   for i=1:nx+1
      @views Ul, Ur  = U[:,i-1], U[:,i]
      f    = num_flux(eq, lam, Ul, Ur, xf[i])
      @views res[:, i-1] += f/ dx0[i-1]
      @views res[:, i]   -= f/ dx0[i]
   end
end

function solve(equation, problem, scheme, param)
   grid = make_grid(problem, param)
   compute_exact_soln! = equation["compute_exact_soln!"]
   nvar = problem["nvar"]
   Tf = problem["final_time"]
   nx = grid.nx
   dx = grid.dx
   xf = grid.xf
   # Allocating variables
   U   = gArray(nvar, nx)
   Ue  = zeros(nvar, nx)
   res = gArray(nvar, nx) # dU/dt + res(U) = 0
   Ua  = U # ua is just Ua for this first order method,
           # storing for clarity
   set_initial_value!(grid, U, problem)
   it, t = 0, 0.0
   figure(figsize=(15,5))
   while t < Tf
      lam, dt = compute_lam_dt(equation, grid, Ua)
      adjust_time_step(problem, param, dt, t)
      compute_exact_soln!(equation["eq"], grid, t, problem, nvar, Ue)
      update_ghost!(grid, U, Ue)                            # Fills ghost cells
      compute_residual!(equation, grid, lam, U, scheme, res)
      @. U -= dt*res
      plot_solution(grid, U, Ue, t, it, param)
      t += dt; it += 1
   end
end

export Problem
export Parameters
export Scheme
export solve

end