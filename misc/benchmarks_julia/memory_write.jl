@inbounds @inline function flux!(F, U)
   γ = 1.4
   u = U[2] / U[1]    # velocity
   p = (γ - 1.0) * (U[3] - 0.5 * U[1] * u^2) # pressure
   F[1], F[2], F[3] = [U[2], p + U[1] * u^2, (U[3]+p) * u] # flux
   return nothing
end

@inline function compute_flux(F, U, nx)
   for i in 1:nx
      @views flux!(F[:,i], U[:,i])
   end
end

function setup_arrays(nvar, nx)
   U = ones(nvar, nx)
   F = Array{Float64}(undef, nvar, nx)
   return U, F
end

function benchmark_in_place(nvar, nx)
   U, F = setup_arrays(nvar, nx)
   @time compute_flux(F, U, nx)
   return F
end
F1 = benchmark_in_place(3, 100) # compile
F1 = benchmark_in_place(3, 100) # actual measurement

