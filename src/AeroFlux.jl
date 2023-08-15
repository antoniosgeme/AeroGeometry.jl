module AeroFlux

using Reexport
include("./Airfoils.jl")
@reexport using .Airfoils

include("./Wings.jl")
@reexport using .Wings

include("./Singularities.jl")
@reexport using .Singularities

end # module AeroFlux

