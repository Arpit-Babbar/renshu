ny = 10
nd = 5
t = 0.0
x = 0.0
dt = 0.5
Ub, Fb = zeros(ny), zeros(ny)
y = 100
@inbounds Threads.@threads for j=1:ny
    for k=1:nd
       y = j+k
       ub, fb = 0.0, 0.0
       for l=1:nd
          tq = t + l * dt
          bvalue = (tq+y)^2
          ub += bvalue * l
          fb += (bvalue-1)*l
       end
       Ub[j] = ub
       Fb[j] = fb
   end
end

println("Ub = ", Ub)
println("y = ", y)
println("Fb = ", Fb)
