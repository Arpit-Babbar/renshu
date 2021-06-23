using Plots
# Source- https://stackoverflow.com/questions/55794068/animating-subplots-using-plots-jl-efficiently

function main()
   x  = LinRange(0.0,1.0,10)
   dt = 0.1
   Tf = 10.0
   t  = 0.0
   u  = sin.(2.0*π*x.-t)
   anim=Animation()
   while t<Tf
      t += dt
      u  = sin.(2.0*π*(x.-t))
      plot(x,u)
      frame(anim)
   end
   plot(p) # final plot
   gif(anim, "anim_fps15.gif", fps = 5)
end
main()
