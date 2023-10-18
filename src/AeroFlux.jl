module AeroFlux

include("./Airfoils.jl")
export Airfoil, plotme, get_upper_coordinates, get_lower_coordinates, get_area,
       get_centroid, repanel!, write_file, get_local_camber, get_local_thickness,
       get_max_camber, get_max_thickness, get_LE_index, get_TE_thickness,
       get_TE_angle, add_control_surface!, repanel, blend_airfoils

include("./Wings.jl")
export Wing, WingXSec

#include("./Singularities.jl")

#@reexport using .Singularities
#export LineVortex

end # module AeroFlux

