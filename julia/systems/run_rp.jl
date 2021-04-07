# Solver for Riemann problem of Linear Advection
# User chooses Ul, Ur

push!(LOAD_PATH,".")
using main
using VTK_OUT
using EqLinAdv

nx = 50 # number of cells

# Plan - To not have to worry about periodicity, we can temporarily
# specify the boundary points to be the Dirichlet values.
# this will work until the wave moves to the boundary, so we will
# test till that short time.

xmin, xmax   = 0.0, 1.0 # domain
nvar         = 3        # number of variables
fprime(U,x) = Diagonal([1.0,2.0,3.0]) # The matrix given by F'(U)
# initial condition
Ul, Ur       = [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]
initial_condition(x) = (x <= 0) ? Ul : Ur
num_flux     = lax_friedrich


# -----------------------------------
eq = LinAdv(fprime)