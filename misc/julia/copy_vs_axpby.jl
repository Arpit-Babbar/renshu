using BenchmarkTools
using LinearAlgebra

n1,n2,nx,ny = 10, 10, 100, 100
u0, v   = rand(n1,n2,nx,ny), rand(n1,n2,nx,ny)
u1       = rand(n1,n2,nx,ny)

function extra_copy!(u0, v, u1)
   u1 .= axpby!(1.0, u0, -0.5, v)
   return nothing
end

function extra_steps_copyto!(u0, v, u1)
   copyto!(u0,u1)
   axpy!(-0.5, v, u1)
   return nothing
end

function extra_steps_assign!(u0, v, u1)
   u1  .= u0
   axpy!(-0.5, v, u1)
   return nothing
end

function extra_steps_lowlevel!(u0, v, u1)
   axpby!(1.0, u0, 0.0, u1)
   axpy!(-0.5, v, u1)
   return nothing
end

println("Old version")
@btime extra_copy!(u0, v, u1)
println("New version with copyto!")
@btime extra_steps_copyto!(u0, v, u1)
println("New version with assignment")
@btime extra_steps_copyto!(u0, v, u1)
println("New version with low level")
@btime extra_steps_lowlevel!(u0, v, u1)

println("This might make one think that axpby! is better than copy!")
println("But, axpby! is actually using threads, so it's not a fair comparison.")
println("If we run the same tests with number of threads set to 1, we get")
BLAS.set_num_threads(1)
println("Old version")
@btime extra_copy!(u0, v, u1)
println("New version with copyto!")
@btime extra_steps_copyto!(u0, v, u1)
println("New version with assignment")
@btime extra_steps_copyto!(u0, v, u1)
println("New version with low level")
@btime extra_steps_lowlevel!(u0, v, u1)