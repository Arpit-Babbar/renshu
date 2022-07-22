using Symbolics
using LinearAlgebra

xmin = -3.0
xmax = 0.0
L = xmax-xmin

N = 6
h = L/(N-1.0)

# @variables h

# x = [xmin + h*i for i=0:N-1]
x = [0.0 for i=0:100]

# A = zeros(N-1,N-1)
A = zeros(Num, N-1,N-1)

A[1,1] = -2.0/h^2 + x[1]^2
A[1,2] =  1.0/h^2

for i=2:N-2
   A[i,i-1] = -1.0/h^2
   A[i,i]   = 2.0/h^2 + x[i]^2
   A[i,i+1] = -1.0/h^2
end
A[N-1,N-3] = -1.0/h^2
A[N-1,N-2] = 2.0/h^2 + x[N-1]^2
A[N-1,N-1] = -1.0/h^2

# println(eigvals(A))
println(simplify(det(A)))
x[N-1]
