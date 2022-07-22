# A * B'
using LinearAlgebra
using BenchmarkTools
nvar = 3
nd = 4
A = rand(nd,nd)
B = rand(nvar,nd)
res = zeros(nvar,nd)

function test(C, A, B)
   mul!(C,B,A')
end

@btime test($res, $A, $B)
function mulT!(C, A, B, a, n1, n2)
   for ii=1:n2
      for k=1:n2
         for n=1:n1
            C[n,k] += a*A[k,ii] * B[n,ii]
         end
      end
   end
end

@btime mulT!($res, $A, $B, 2.0, $nvar, $nd)
