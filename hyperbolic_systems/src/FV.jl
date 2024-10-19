module FV

using ..Grid
using LinearAlgebra
using OffsetArrays
using TickTock
using Plots
using UnPack
using StaticArrays
using Cthulhu
# plotly(size = (750, 565)) # use gr as the plot backend, for its high performance
#-------------------------------------------------------------------------------
# Create a dictionary of problem description
#-------------------------------------------------------------------------------
function Problem(domain::Tuple{Float64,Float64},
                 nvar::Int64,
                 initial_value::Function,
                 boundary_value::Function,
                 boundary_condition::Any,
                 final_time::Float64)
   return (; domain, nvar, initial_value, boundary_value,
             boundary_condition, final_time)
end

#-------------------------------------------------------------------------------
# Create a dictionary of parameters
#-------------------------------------------------------------------------------
function Parameters(grid_size::Int64,
                    cfl::Float64,
                    Ccfl::Float64,
                    save_time_interval::Float64)
   @assert (cfl >= 0.0) "cfl must be >= 0.0"
   @assert (save_time_interval >= 0.0) "save_time_interval must be >= 0.0"
   return (;grid_size, cfl, Ccfl, save_time_interval)
end

#-------------------------------------------------------------------------------
# Create a dictionary of scheme description
#-------------------------------------------------------------------------------
function Scheme(equation, numflux::String)
   @unpack numfluxes = equation
   (;numflux = numfluxes[numflux], numflux_ind = numflux)
end

#-------------------------------------------------------------------------------
# Convert array to string
#-------------------------------------------------------------------------------
function array2string(arr)
   arr_string = "["
   n = size(arr)[1]
   for i=1:n-1
      arr_string = arr_string * string(arr[i]) * ","
   end
   arr_string = arr_string * string(arr[end]) * "]"
end

#-------------------------------------------------------------------------------
# Create ghosted arrays with one layer of cells on both sides
#-------------------------------------------------------------------------------
gArray(nvar, nx) = OffsetArray(zeros(nvar, nx+2),
                              OffsetArrays.Origin(1, 0))


#-------------------------------------------------------------------------------
# Adjust dt to reach final time or the next time when solution has to be saved
#-------------------------------------------------------------------------------
# function adjust_time_step_cuda(problem, param, dt, t)
#    # Adjust to reach final time exactly
#    @unpack final_time = problem
#    @unpack save_time_interval = param
#    CUDA.@allowscalar if t + dt[1] > final_time
#       dt[1] = final_time - t
#       return nothing
#    end

#    # Adjust to reach next solution saving time
#    CUDA.@allowscalar if save_time_interval > 0.0
#       next_save_time = ceil(t/save_time_interval) * save_time_interval
#       # If t is not a plotting time, we check if the next time
#       # would step over the plotting time to adjust dt
#       if abs(t-next_save_time) > 1e-10 && t + dt[1] - next_save_time > 1e-10
#          dt[1] = next_save_time - t
#          return nothing
#       end
#    end
#    return nothing
# end

function adjust_time_step(problem, param, dt, t)
   # Adjust to reach final time exactly
   @unpack final_time = problem
   @unpack save_time_interval = param
   if t + dt[1] > final_time
      dt[1] = final_time - t
      return nothing
   end

   # Adjust to reach next solution saving time
   if save_time_interval > 0.0
      next_save_time = ceil(t/save_time_interval) * save_time_interval
      # If t is not a plotting time, we check if the next time
      # would step over the plotting time to adjust dt
      if abs(t-next_save_time) > 1e-10 && t + dt[1] - next_save_time > 1e-10
         dt[1] = next_save_time - t
         return nothing
      end
   end
   return nothing
end

function compute_dt(equation, scheme, param, grid, Ua, dt, lam)
   @unpack fprime, eq = equation
   @unpack Ccfl = param
   nx, dx = grid.nx, grid.dx
   xc     = grid.xc
   lam_loc = 0.0
   dt_loc  = 1.0
   for i=1:nx
      ua   = Ua[:, i]
      lam0 = maximum(abs.(eigvals(equation.fprime(xc[i], ua, eq)))) # CHECK xc[i]!
      lam_loc  = max(lam_loc, lam0)
      dt_loc   = min(dt_loc, dx[i]/lam0)
   end
   dt[1]   = dt_loc
   lam[1] = lam_loc[1]
   dt .= Ccfl*dt
end

function compute_dt_cuda(eq, grid, Ua, dt, lam)

   # CUDA.sync_warp()
   return nothing
end

function set_initial_value!(grid, U, initial_value)
   nx = grid.nx
   xc = grid.xc
   for i=1:nx
      @views initial_value(U[:,i], xc[i])
   end
end

function set_initial_value_cuda!(grid, U, initial_value)
   nx = grid.nx
   xc = grid.xc
   idx = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
   # idx = 1
   i = idx+1
   if i > nx+1
      return nothing
   end
   u = @view U[:,i]
   @views initial_value(u, xc[i-1])
   return nothing
end

function update_ghost!(grid, U, problem)
   xmin, xmax = grid.domain
   nx = grid.nx
   initial_value = problem.initial_value
   if problem.boundary_condition == "Dirichlet"
      @views U[:, 0]   .= U[:, 1] # only for short time
      @views U[:,nx+1] .= U[:,nx]
   else
      @views U[:, 0]    .= U[:, nx]
      @views U[:, nx+1] .= U[:,1]
   end
   return nothing
end

function my_fill!(uf, init)
   idx = CUDA.threadIdx().x
   uf[idx] = init
   return nothing
end

@inline function get_node_vars(u,indices...)
   SVector(ntuple(@inline(v -> u[v, indices...]), 3))
end

# function compute_residual_cuda!(equation, grid, lam, U, scheme, res,
#                            dx0, Uf, promoter)
#    # TODO - Do we really want one thread to do just one iteration?
#    nx = grid.nx
#    xf = grid.xf
#    dx = grid.dx
#    @unpack eq = equation
#    @unpack numflux = scheme
#    dx0[2:nx+1] .= dx # std array
#    # res[:,:] .= 0.0 # Shouldn't we be able to avoid this?
#    fill!(res, zero(eltype(res)))
#                    # Something like this?
#                    #      @views res[:, i-1] += f/ dx0[i-1]
#                    #      @views res[:, i]   = f/ dx0[i]
#    # loop over faces
#    @cuda threads=3 my_fill!(Uf, zero(eltype(Uf)))
#    @unpack eq = equation
#    γ = eq.γ
#    function cuda_res1(γ, lam, U, xf, Uf, res, dx0, nx)
#       Uf[1] = 0.0
#       Uf[2] = 0.0
#       Uf[3] = 0.0
#       # Ul = get_node_vars(U, i-1)
#       idx = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
#       i = idx+1
#       if i > nx+2
#          return nothing
#       end
#       Ul = SVector{3}(U[1,i-1], U[2,i-1], U[3,i-1])
#       Ur = SVector{3}(U[1,i], U[2,i], U[3,i])
#       Uf_ = numflux(γ, lam[1], Ul, Ur, xf[i-1], Uf)
#       for n in 1:3
#          res[n, i-1] += Uf_[n] / dx0[i-1]
#       end
#       return nothing
#    end

#    function cuda_res2(γ, lam, U, xf, Uf, res, dx0, nx)
#       Uf[1] = 0.0
#       Uf[2] = 0.0
#       Uf[3] = 0.0
#       idx = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
#       # idx = 1
#       i = idx+1
#       if i > nx+2
#          return nothing
#       end
#       Ul = SVector{3}(U[1,i-1], U[2,i-1], U[3,i-1])
#       Ur = SVector{3}(U[1,i], U[2,i], U[3,i])
#       Uf_ = numflux(γ, lam[1], Ul, Ur, xf[i-1], Uf)
#       for n in 1:3
#          res[n, i]   -= Uf_[n] / dx0[i]
#       end
#       return nothing
#    end

#    nthreads = min(400, nx+1)
#    nblocks  = ceil(Int64, ( (nx+1)/nthreads) )

#    @cuda threads=nthreads blocks=nblocks cuda_res1(γ, lam, U, xf, Uf, res, dx0, nx) # std array
#    @cuda threads=nthreads blocks=nblocks cuda_res2(γ, lam, U, xf, Uf, res, dx0, nx) # std array
#    return nothing
# end

function compute_error(grid, U, t, equation, problem)
   @unpack boundary_value = problem
   error_l2 = 0.0
   error_l1 = 0.0
   error_linf = 0.0
   nx, dx = grid.nx, grid.dx
   xc     = grid.xc
   lam_loc = 0.0
   dt_loc  = 1.0
   for i=1:nx
      u_   = U[1, i]
      u_exact = boundary_value(xc[i], t)
      error = abs(u_ - u_exact[1])
      error_l1 += error   * dx[i]
      error_l2 += error^2 * dx[i]
      error_linf = max(error_linf, error)
   end
   return error_l1, error_l2, error_linf
end

function compute_residual!(equation, grid, lam, U, scheme, res,
                           dx0, Uf, promoter)
   nx = grid.nx
   xf = grid.xf
   dx = grid.dx
   eq = equation.eq
   numflux = scheme.numflux
   dx0 =  OffsetArray(zeros(nx+2), OffsetArrays.Origin(0))
   dx0[1:nx] .= dx
   dx0[0] = dx0[nx+1] = 0.0 # redundant values
   res[:,:] .= 0.0 # Shouldn't we be able to avoid this?
                   # Something like this?
                   #      @views res[:, i-1] += f/ dx0[i-1]
                   #      @views res[:, i]   = f/ dx0[i]
   # loop over faces
   Uf = zeros(3)
   for i=1:nx+1
      @views Ul, Ur  = U[:,i-1], U[:,i]
      numflux(equation.eq, lam, Ul, Ur, xf[i], Uf)
      @views res[:, i-1] += Uf/ dx0[i-1]
      @views res[:, i]   -= Uf/ dx0[i]
   end
end

# function solve_cuda(equation, problem, scheme, param, promoter)
#    tick()
#    grid = make_grid(problem, param, promoter)
#    @unpack nvar = problem
#    @unpack final_time = problem
#    Tf = final_time
#    nx = grid.nx
#    # Allocating variables

#    U   = CuArray(zeros(nvar, nx+2)) # std array
#    res = CuArray(zeros(nvar, nx+2)) # dU/dt + res(U) = 0
#    dx0 =  CuArray(zeros(nx + 2))
#    Uf = CuArray(zeros(3))
#    dt, lam = CuArray([1e20]), CuArray([0.0])

#    Ua  = U
#    @unpack initial_value = problem
#    set_initial_value!(grid, U, initial_value)
#    it, t = 0, 0.0
#    while t < Tf
#       compute_dt(equation, scheme, param, grid, Ua, dt, lam)
#       adjust_time_step(problem, param, dt, t)
#       update_ghost!(grid, U, problem)
#       compute_residual!(equation, grid, lam, U, scheme, res,
#                         dx0, Uf, promoter)
#       @. U -= dt*res
#       CUDA.@allowscalar t += dt[1]; it += 1
#    end
#    tock()
#    sol = (grid, U)
#    return sol
# end

function solve(equation, problem, scheme, param, promoter, plotters)
   tick()
   grid = make_grid(problem, param, promoter)
   initialize_plot = plotters["initialize_plot"]
   update_plot!    = plotters["update_plot!"]
   final_plot      = plotters["final_plot"]
   @unpack nvar = problem
   @unpack final_time = problem
   Tf = final_time
   nx = grid.nx
   # Allocating variables

   U   = gArray(nvar, nx)
   Ue  = zeros(nvar, nx)
   res = gArray(nvar, nx) # dU/dt + res(U) = 0
   Ua  = U # ua is just Ua for this first order method,
           # storing for clarity
   dx0 =  zeros(nx + 2)
   Uf = zeros(3)
   dt, lam = [1e20], [0.0]

   @unpack initial_value = problem
   set_initial_value!(grid, U, initial_value)
   it, t = 0, 0.0
   plt_data = initialize_plot(grid, problem, equation.eq, scheme, U)
   while t < Tf
      l1, l2, linf = compute_error(grid, U, t, equation, problem)
      @show l1, l2, linf
      compute_dt(equation, scheme, param, grid, Ua, dt, lam)
      adjust_time_step(problem, param, dt, t)
      update_ghost!(grid, U, problem)
      compute_residual!(equation, grid, lam, U, scheme, res,
                        dx0, Uf, promoter)
      @. U -= dt*res
      t += dt[1]; it += 1
      @show t, dt
      update_plot!(grid, problem, equation.eq, scheme, U, t, it, param, plt_data)
   end
   l1, l2, linf = compute_error(grid, U, t, equation, problem)

   final_plot(plt_data, equation.eq)
   p, anim = plt_data

   tock()
   sol = (;p, anim, grid, U, l1, l2, linf)
   return sol
end

export Problem
export Parameters
export Scheme
export empty_func
export array2string
export solve

end