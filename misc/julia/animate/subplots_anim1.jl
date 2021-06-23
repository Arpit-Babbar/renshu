using Plots
# Source- https://stackoverflow.com/questions/55794068/animating-subplots-using-plots-jl-efficiently

# functions to test plotting
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

function generate_plot(x, u)
   plot(x,u)
end

function update_plot(x, u)
   plot(x,u)
end

function main()
   x  = LinRange(0.0,1.0,10)
   dt = 0.1
   Tf = 10.0
   t  = 0.0
   l = @layout[ a b  c]
   p1 = plot(x,sine.(x.-t),label="sine")
   p2 = plot(x,flat_hat.(x.-t),label="flat_hat")
   p3 = plot(x,linear.(x.-t),label="linear",ylim=(-1,1))
   p = plot(p1, p2, p3, layout = l)
   anim = @animate while t<Tf
      t += dt
      p[1][1][:y] = sin.(x.-t)
      p[2][1][:y] = flat_hat.(x.-t)
      p[3][1][:y] = linear.(x.-t)
   end
   plot(p) # final plot
   gif(anim, "anim_fps15.gif", fps = 5)
end
main()
