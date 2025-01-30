module AeroGeometry
using LinearAlgebra
using Statistics

include("./Tools.jl")


include("./Airfoils.jl")
export Airfoil, coordinates, area,
       centroid, repanel!, write_file, local_camber, local_thickness,
       max_camber, max_thickness, leading_edge_index, trailing_edge_thickness,
       trailing_edge_angle, deflect_control_surface!,deflect_control_surface , repanel, blend_airfoils,
       surface_coordinates, quarter_chord

include("./Wings.jl")
export Wing, WingXSec, ControlSurface, translate! 

include("./Fuselages.jl")
export Fuselage, FuselageXSec, coordinates, add_xsec!

include("./Airplanes.jl")
export Airplane

include("./Visualizations.jl")

include("./Library.jl")
export Cessna152

end