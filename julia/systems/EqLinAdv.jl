module EqLinAdv

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

function upwind(eq, lam, Ul, Ur, x) # Numerical flux of face at x
   Fl, Fr = flux(x, Ul, eq), flux(x, Ur, eq)
   value  = Fl
   return value
end

export LinAdv
export flux
export lax_friedrich
export upwind

end