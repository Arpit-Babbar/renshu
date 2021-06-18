using FLoops
n = 10

a = LinRange(1,n,n)
c = 0.0
d = 0.0
println("Number of threads = ", Threads.nthreads())
@floop ThreadedEx(basesize=4) for i in 2:n-1
   d1 = a[i]
   @reduce(c += d1)
   d2 = a[i-1]
   @reduce(c += d2)
   d3 = a[i+1]
   @reduce(c += d3)
   println(Threads.threadid())
end
c
println(c)
