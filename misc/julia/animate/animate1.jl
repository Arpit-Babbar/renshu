using Plots

function generate_plot(x, u)
   plot(x,u)
end

function update_plot(x,u)
   plot(x,u)
end

x  = LinRange(0.0,1.0,10)
dt = 0.1
Tf = 10.0
t  = 0.0
u  = sin.(2.0*π*x.-t)
p = plot(x,u)
anim=Animation()
anim = @animate while t<Tf
   t += dt
   u  = sin.(2.0*π*(x.-t))
   p = plot(x,u)
   frame(anim)
end
gif(anim, "anim_fps15.gif", fps = 5)
wait(10)
# plot(x,u)