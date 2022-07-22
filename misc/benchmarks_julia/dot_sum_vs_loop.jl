# A * B'
using LinearAlgebra
using BenchmarkTools
nvar = 3
# nd = 4
A = rand(nvar)
B = rand(nvar)
res = zeros(nvar)

function sum_loop!(A,B,res, nvar)
   for n=1:nvar
      res[n] = 3.0*A[n]+4.0*B[n] + (A[n]-B[n])/0.9
   end
end

function sum_dot!(A, B, res)
   @. res = 3.0*A+4.0*B + (A-B)/0.9
end

@btime sum_loop!($A,$B,$res,$nvar)

@btime sum_dot!($A,$B,$res)
