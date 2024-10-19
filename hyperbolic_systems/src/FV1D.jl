module FV1D

function get_node_vars(U, eq, indices)
    SVector(ntuple(v -> U[v, indices], 1))
 end

include("Grid.jl")
include("FV.jl")
include("EqEuler.jl")
include("EqIsentropicEuler.jl")
include("EqLinAdv.jl")

end # module
