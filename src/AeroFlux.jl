module AeroFlux

using Reexport

include("./airfoil.jl")
@reexport using .airfoil

foil = Airfoil("naca6412")
add_control_surface!(foil,deflection=10)

end # module AeroFlux

