using Plots

function sine(x)
   return sin(2.0*π*x)
end

function flat_hat(x)
   x0 = x % 1.0
   if 0.25<=abs(x0)<=0.75
      return 1.0
   else
      return 0.0
   end
end

function linear(x)
   x0 = x % 1.0
   return x0
end

function main()
   x  = LinRange(0.0,1.0,10)
   anim=Animation()
   l = @layout[ a b c]
   p1 = plot(x,sine.(x.-t),label="sine")
   p2 = plot(x,flat_hat.(x.-t),label="flat_hat")
   p3 = plot(x,linear.(x.-t),label="linear",ylim=(-1,1))
   p = plot(p1, p2, p3, layout = l)
end

main()
