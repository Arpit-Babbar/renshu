module EqLinAdv

struct LinAdv
    fprime::Function
end

flux(x, U, eq::LinAdv) = eq.fprime(U,x) * U
                        # matrix multiplication

function lax_friedrich(eq, lam, Ul, Ur, x) # Numerical flux of face at x
    Fl, Fr = eq.flux(x, Ul, eq), eq.flux(x, Ur, eq)
    value  = 0.5*(Fl+Fr) - 0.5 * lam * (Ur - Ul)
    return value
end

end