module AeroFlux

using Reexport
#using Plots

include("./Airfoils.jl")
@reexport using .Airfoils

#foil = Airfoil("naca6412")
#fig = plotme(foil)
#repanel!(foil,100)
#scatter!(fig,foil.coordinates[:,1],foil.coordinates[:,2])
#add_control_surface!(foil,deflection=10)

end # module AeroFlux

