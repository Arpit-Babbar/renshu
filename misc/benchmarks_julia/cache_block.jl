using StaticArrays


abstract type AbstractEquations{NDIMS, NVARS} end

struct Euler <: AbstractEquations{1, 3}
end

@inline nvariables(::AbstractEquations{NDIMS,NVARS}) where {NDIMS, NVARS} = NVARS

@inline function get_node_vars(u, eq, indices...)
   SVector(ntuple(@inline(v -> u[v, indices...]), Val(nvariables(eq))))
end

@inline function flux(u)
   gamma = 1.4
   rho, rho_v1, rho_e = u
   v1 = rho_v1 / rho
   p  = (gamma - 1) * (rho_e - 0.5 * rho_v1 * v1)
   f1 = rho_v1
   f2 = rho_v1 * v1 + p
   f3 = (rho_e + p) * v1
   return SVector(f1, f2, f3)
end

@inline function set_node_vars!(u, u_node, eq, indices...)
   for v in Base.OneTo(nvariables(eq))
      u[v, indices...] = u_node[v]
   end
   return nothing
end

function compute_flux_trixi(F, U, nvar, nx)
   for i in 1:nx
      u_node = get_node_vars(U, nvar, i)
      f_node = flux(u_node)
      set_node_vars!(F, f_node, nvar, i)
   end
end

function setup_arrays(eq::Euler, nx)
   nvar = nvariables(eq)
   U = ones(nvar, nx)
   F = Array{Float64}(undef, nvar, nx)
   return U, F
end

function benchmark_trixi(nx)
   eq = Euler()
   U, F = setup_arrays(eq, nx)
   @time compute_flux_trixi(F, U, eq, nx)
   return F
end
benchmark_trixi(100); # compile
benchmark_trixi(100); # actual measurement

