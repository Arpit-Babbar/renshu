using BenchmarkTools

function foo1(a)
   a = a.^2
   return nothing
end

function foo2(a)
   @. a = a^2
   return nothing
end

function foo3(a)
   a .= a.^2
   return nothing
end

n = 10000
a1 = rand(n)
@btime foo1(a1)
a2 = rand(n)
@btime foo2(a2)
a3 = rand(n)
@btime foo3(a3)

# The output we get is
# 3.117 μs (2 allocations: 78.20 KiB)
# 1.010 μs (0 allocations: 0 bytes)
# 1.080 μs (0 allocations: 0 bytes)
# So, the best methods are
# @. a = a^2
# or
# a. = a.^2
