module Grid

using Printf

struct CartesianGrid
   domain::Tuple{Float64,Float64}  # xmin, xmax
   nx::Int64                 # nx - number of points
   xc::Array{Float64,1}      # cell centers
   xf::Array{Float64,1}      # cell faces
   dx::Array{Float64,1}      # cell sizes
end

# Uniform Cartesian grid - Change this function to run solver for different
#                          grids
function make_grid(problem, param)
   xmin, xmax = problem["domain"]
   nx         = param["grid_size"]
   println("Making uniform grid of interval [", xmin, ", ", xmax,"]")
   dx1 = (xmax - xmin)/nx
   xc = LinRange(xmin+0.5*dx1, xmax-0.5*dx1, nx)
   @printf("   Grid of with number of points = %d \n", nx)
   @printf("   xmin,xmax                     = %e, %e\n", xmin, xmax)
   @printf("   dx                            = %e\n", dx1)
   dx = dx1 .* ones(nx)
   xf = LinRange(xmin, xmax, nx+1)
   return CartesianGrid((xmin, xmax), nx, xc, xf, dx)
end

export make_grid

end