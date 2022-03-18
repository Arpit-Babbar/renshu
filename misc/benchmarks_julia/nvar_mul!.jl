using LoopVectorization
using BenchmarkTools
using LinearAlgebra
using StaticArrays

# @inline function mult!(C, A, B, nvar, nd)
#    fill!(C,zero(eltype(C)))
#    # Default sets thread = false
#    @turbo for j=1:nd, k=1:nd, n=1:nvar
#       for i=1:nd
#          C[n,i,j] += A[i,k]*B[n,k,j]
#       end
#    end
# end

@inline function mult!(C, A, B)
   fill!(C, zero(eltype(C)))
   @turbo for n ∈ axes(C,1), i ∈ axes(A,2), j ∈ axes(B,3)
      for k ∈ axes(A,2)
        C[n,i,j] += A[i,k] * B[n,k,j]
      end
   end
end

@inline function mygemmavx!(C, A, B)
   fill!(C, zero(eltype(C)))
   @turbo for n ∈ axes(A,1), m ∈ axes(B,2)
      #  Cmn = zero(eltype(C))
      for k ∈ axes(A,2)
        C[m,n] += A[m,k] * B[k,n]
      end
      #  C[m,n] = Cmn
   end
end

nd = 4
nvar = 1
# C = zeros(nvar, nd, nd)
# B = rand(nvar, nd, nd)
# A = rand(nd, nd)

C = zeros(nvar, nd, nd)
B = @SArray rand(nvar, nd, nd)
A = @SArray rand(nvar,nd, nd)

C_ = zeros(nd, nd)
B_ = @SArray rand(nd, nd)
A_ = @SArray rand(nd, nd)

@btime mult!($C, $A, $B);

# @btime mul!($C_, $A_, $B_);
@btime mygemmavx!($C_, $A_, $B_);
