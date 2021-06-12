using LinearAlgebra
using BenchmarkTools

function loops(u1, Vl, Vr, wg, nd)
   ul, ur, ub, ut = 0.0, 0.0, 0.0, 0.0
   for jj=1:nd, ii=1:nd
      ul += u1[jj, ii, 1, 1] * Vl[jj] * wg[ii] # transpose(u1) * Vl . wg
      ur += u1[jj, ii, 1, 1] * Vr[jj] * wg[ii] # transpose(u1) * Vr . wg
      ub += u1[ii, jj, 1, 1] * Vl[jj] * wg[ii] # u1 * Vl * wg
      ut += u1[ii, jj, 1, 1] * Vr[jj] * wg[ii] # u1 * Vl * wg
   end
   return ul, ur, ub, ut
end

function low_level(u1, Vl, Vr, wg, nd, Ub)
   @views U, UT = u1[:,:,1,1], transpose(u1[:,:,1,1])
   mul!(Ub, UT, Vl)
   ul = BLAS.dot(nd, Ub, 1, wg, 1)
   mul!(Ub, UT, Vr)
   ur = BLAS.dot(nd, Ub, 1, wg, 1)
   mul!(Ub, U,  Vl)
   ub = BLAS.dot(nd, Ub, 1, wg, 1)
   mul!(Ub, U,  Vr)
   ut = BLAS.dot(nd, Ub, 1, wg, 1)
   return ul, ur, ub, ut
end

nx = ny = 100
nd = 5
Vl, Vr, wg = rand(nd), rand(nd), rand(nd)
u1 = rand(nd,nd,nx,ny)

@btime loops(u1, Vl, Vr, wg, nd)
Ub = zeros(nd)
@btime low_level(u1, Vl, Vr, wg, nd, Ub)

# nd = 2
# loops - 32.998 ns (1 allocation: 48 bytes)
# low_level - 187.278 ns (1 allocation: 48 bytes)
# nd = 3
# loops - 44.747 ns (1 allocation: 48 bytes)
# low_level - 189.628 ns (1 allocation: 48 bytes)
# nd = 3
# loops - 45.354 ns (1 allocation: 48 bytes)
# low_level - 195.130 ns (1 allocation: 48 bytes)
# nd = 4
# loops - 62.347 ns (1 allocation: 48 bytes)
# low-level - 214.499 ns (1 allocation: 48 bytes)
# nd = 5
# loops - 87.161 ns (1 allocation: 48 bytes)
# low-level - 246.797 ns (1 allocation: 48 bytes)
# nd = 100
# loops - 24.500 μs (1 allocation: 48 bytes)
# low-level - 6.250 μs (1 allocation: 48 bytes)
