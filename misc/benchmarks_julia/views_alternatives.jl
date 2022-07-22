using BenchmarkTools
using StaticArrays

struct StructWithFunc1
   f::Function
end

struct StructWithFunc2{F <: Function}
   f::F
end

@inbounds @inline function random_math_viewed(d,a,al,ar,b,bl,br,c,cl,cr,nvar)
   for n=1:nvar
      d[n] = 0.5*(al[n]+bl[n]) + 129.0*(c[n]*br[n]) + 300.0*(al[n]*cr[n])
   end
end

struct_with_func1 = StructWithFunc1(random_math_viewed)
struct_with_func2 = StructWithFunc2(random_math_viewed)

@inbounds @inline function random_math_indexed(d,a0,b0,c0, i, nvar)
   for n=1:nvar
      d[n,i] = 0.5*(a0[n,i]+b0[n,i-1]) + 129.0*(c0[n,i]*b0[n,i+1]) + 300.0*(a0[n,i-1]*c0[n,i+1])
   end
end

@inline function get_node_vars(u, i, nvar)
   SVector(ntuple(@inline(v -> u[v,i]), Val(nvar)))
end

@inline function get_node_vars(u, i, nvar)
   SVector(ntuple(@inline(v -> u[v,i]), Val(nvar)))
end

function set_node_vars!(u, u_node, i, nvar)
   for n=1:nvar
      u[n,i] = u_node[n]
   end
end

@inline @inbounds function random_math_returning(a,al,ar,b,bl,br,c,cl,cr,nvar)
   # Assumes nvar=4, which is fine
   d1 = 0.5*(a[1]+bl[1]) + 129.0*(c[1]*br[1]) + 300.0*(al[1]*cr[1])
   d2 = 0.5*(a[2]+bl[2]) + 129.0*(c[2]*br[2]) + 300.0*(al[2]*cr[2])
   d3 = 0.5*(a[3]+bl[3]) + 129.0*(c[3]*br[3]) + 300.0*(al[3]*cr[3])
   d4 = 0.5*(a[4]+bl[4]) + 129.0*(c[4]*br[4]) + 300.0*(al[4]*cr[4])
   return SVector(d1,d2,d3,d4)
end

nx = 1000000
nvar = 4

a,b,c,d = [rand(nvar,nx) for _=1:4]
function func_indexed(a, b, c, d, nx, nvar)
   for i=2:nx-1
      random_math_indexed(d, a, b, c, i, nvar)
   end
end

function func_viewed(a,b,c,d,nx,nvar)
   @inbounds for i=2:nx-1
      a0, al, ar = @views a[:,i], a[:,i-1], a[:,i+1]
      b0, bl, br  = @views b[:,i], b[:,i-1], b[:,i+1]
      c0, cl, cr  = @views c[:,i], c[:,i-1], c[:,i+1]
      d0  = @views d[:,i]
      random_math_viewed(d0,a0,al,ar,b0,bl,br,c0,cl,cr,nvar)
   end
end

function func_viewed_with_struct(a,b,c,d,nx,nvar,struct_with_func)
   f = struct_with_func.f
   @inbounds for i=2:nx-1
      a0, al, ar = @views a[:,i], a[:,i-1], a[:,i+1]
      b0, bl, br  = @views b[:,i], b[:,i-1], b[:,i+1]
      c0, cl, cr  = @views c[:,i], c[:,i-1], c[:,i+1]
      d0  = @views d[:,i]
      f(d0,a0,al,ar,b0,bl,br,c0,cl,cr,nvar)
   end
end

function func_copied(a,b,c,d,nx,nvar)
   for i=2:nx-1
      a0, al, ar =  a[:,i], a[:,i-1], a[:,i+1]
      b0, bl, br  =  b[:,i], b[:,i-1], b[:,i+1]
      c0, cl, cr  =  c[:,i], c[:,i-1], c[:,i+1]
      d0  =  d[:,i]
      random_math_viewed(d0,a0,al,ar,b0,bl,br,c0,cl,cr,nvar)
   end
end

function get_node_vars(u, i, nvar)
   SVector(ntuple(@inline(v -> u[v,i]), Val(nvar)))
end

function set_node_vars!(u, u_node, i, nvar)
   for n=1:nvar
      u[n,i] = u_node[n]
   end
end

function func_set_var(a, b, c, d, nx, nvar)
   for i=2:nx-1
         a0, al, ar = ( get_node_vars(a, i, nvar), get_node_vars(a, i-1, nvar),
                        get_node_vars(a, i+1, nvar) )
         b0, bl, br = ( get_node_vars(b, i, nvar), get_node_vars(b, i-1, nvar),
                        get_node_vars(b, i+1, nvar) )
         c0, cl, cr = ( get_node_vars(c, i, nvar), get_node_vars(c, i-1, nvar),
                        get_node_vars(c, i+1, nvar) )
         d0 = random_math_returning(a,al,ar,b,bl,br,c,cl,cr,nvar)
         set_node_vars!(d, d0, i, nvar)
   end

@btime func_indexed($a,$b,$c,$d,$nx,$nvar)
@btime func_viewed($a,$b,$c,$d,$nx,$nvar)
@btime func_viewed_with_struct($a,$b,$c,$d,$nx,$nvar,struct_with_func1)
@btime func_viewed_with_struct($a,$b,$c,$d,$nx,$nvar,struct_with_func2)
# @btime func_copied($a,$b,$c,$d,$nx,$nvar)
@btime func_set_var($a, $b, $c, $d, $nx, $nvar)

# test_performance()