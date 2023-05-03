module AeroFlux

using Reexport
include("./Airfoils.jl")
@reexport using .Airfoils

end # module AeroFlux

