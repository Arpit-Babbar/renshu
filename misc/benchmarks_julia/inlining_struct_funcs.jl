using BenchmarkTools
using StaticArrays


struct StructWithFunc{F <: Function}
   f::F
end

@inbounds @inline function random_math_viewed(d,a,al,ar,b,bl,br,c,cl,cr,nvar)
   for n=1:nvar
      d[n] = 0.5*(al[n]+bl[n]) + 129.0*(c[n]*br[n]) + 300.0*(al[n]*cr[n])
   end
end

struct_with_func1 = StructWithFunc(random_math_viewed)

struct StructWithFunc2{F1,F2 <: Function}
   f1::F1
   f2::F2
end

@inbounds @inline function random_math_viewed1(d,a0,al,ar,b0,bl,br,c0,cl,cr,nvar,
                                               struct_with_func)
   for n=1:nvar
      a, b, c = al[n]+bl[n], c0[n]*br[n], al[n]*cr[n]
      d[n] = struct_with_func.f2(a, b, c)
   end
end

@inbounds @inline function random_math_viewed2(a,b,c)
   return 0.5*a + 129.0*b + 300.0*c
end

struct_with_func2 = StructWithFunc2(random_math_viewed1,random_math_viewed2)

nx = 1000000
nvar = 4

a,b,c,d = [rand(nvar,nx) for _=1:4]
function func_indexed(a, b, c, d, nx, nvar)
   for i=2:nx-1
      random_math_indexed(d, a, b, c, i, nvar)
   end
end

function func_viewed_with_struct(a,b,c,d,nx,nvar,f)
   @inbounds for i=2:nx-1
      a0, al, ar = @views a[:,i], a[:,i-1], a[:,i+1]
      b0, bl, br  = @views b[:,i], b[:,i-1], b[:,i+1]
      c0, cl, cr  = @views c[:,i], c[:,i-1], c[:,i+1]
      d0  = @views d[:,i]
      f(d0,a0,al,ar,b0,bl,br,c0,cl,cr,nvar)
   end
end

# function func_viewed_broken(a,b,c,d,nx,nvar, f1)
#    @inbounds for i=2:nx-1
#       a0, al, ar = @views a[:,i], a[:,i-1], a[:,i+1]
#       b0, bl, br  = @views b[:,i], b[:,i-1], b[:,i+1]
#       c0, cl, cr  = @views c[:,i], c[:,i-1], c[:,i+1]
#       d0  = @views d[:,i]
#       random_math_viewed1(d0,a0,al,ar,b0,bl,br,c0,cl,cr,nvar)
#    end
# end

function func_viewed_broken_with_struct(a,b,c,d,nx,nvar,struct_with_func2)
   f1 = struct_with_func2.f1
   @inbounds for i=2:nx-1
      a0, al, ar = @views a[:,i], a[:,i-1], a[:,i+1]
      b0, bl, br  = @views b[:,i], b[:,i-1], b[:,i+1]
      c0, cl, cr  = @views c[:,i], c[:,i-1], c[:,i+1]
      d0  = @views d[:,i]
      f1(d0, a0, al, ar, b0, bl, br, c0, cl, cr, nvar, struct_with_func2)
   end
end

f = struct_with_func1.f
@btime func_viewed_with_struct($a,$b,$c,$d,$nx,$nvar,$f)
f1 = struct_with_func2.f1
@btime func_viewed_broken_with_struct($a,$b,$c,$d,$nx,$nvar,$struct_with_func2)

# test_performance()
