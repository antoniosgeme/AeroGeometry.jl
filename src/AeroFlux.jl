module AeroFlux

using Reexport
include("./Airfoils.jl")
@reexport using .Airfoils

include("./Singularities.jl")
@reexport using .Singularities

end # module AeroFlux

