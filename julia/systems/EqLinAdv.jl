module EqLinAdv

using LinearAlgebra

struct LinAdv
   fprime::Function
end

flux(x, U, eq::LinAdv) = eq.fprime(U,x) * U
                        # matrix multiplication

function lax_friedrich(eq, lam, Ul, Ur, x) # Numerical flux of face at x
   Fl, Fr = flux(x, Ul, eq), flux(x, Ur, eq)
   value  = 0.5*(Fl+Fr) - 0.5 * lam * (Ur - Ul)
   return value
end

# Looks expensive for a non-constant speed. Can it be fixed?
function upwind(eq, lam, Ul, Ur, x) # Numerical flux of face at x
   eigen_decomp   = eigen(eq.fprime(0.5 * (Ul+Ur), x)) # Roe's scheme had (f(Ul)-f(Ur))/(Ul-Ur)
                                 # appearing hard to extend to high dimensions because Ul, Ur
                                 # are now vectors and can't be put in denominators.
   # for linear advection, it is working because fprime is a constant matrix.
   evalues, evecs = eigen_decomp.values, eigen_decomp.vectors
   lp, lm     = max.(evalues, 0.0), min.(evalues, 0.0)
   lamp, lamm = diagm(lp), diagm(lm)
   Fp = evecs * lamp * inv(evecs) * Ul
   Fm = evecs * lamm * inv(evecs) * Ur
   value  = Fp + Fm
   return value
end

function compute_exact_soln!(eq, grid, t, problem, nvar, Ue)
   nx = grid.nx
   xc = grid.xc
   fp = eq.fprime(Ue, 1.0) # 1.0 is dummy
   initial_value = problem["initial_value"]
   eigen_decomp = eigen(fp)
   Lam, eigen_vecs = eigen_decomp.values, eigen_decomp.vectors
   for j=1:nx
      for i=1:nvar
         # Too complicated, fix.
         Ue[i,j] = (inv(eigen_vecs) * initial_value(xc[j] - Lam[i] * t))[i]
      end
   end
   for j=1:nx
      Ue[:,j] = eigen_vecs * Ue[:,j]
   end
   return nothing
end

get_equation(fprime) = Dict("eq"                  => LinAdv(fprime),
                            "flux"                => flux,
                            "speed"               => fprime,
                            "compute_exact_soln!" => compute_exact_soln!,
                            "name"                => "2d Linear advection equation")

export lax_friedrich
export upwind
export get_equation

end