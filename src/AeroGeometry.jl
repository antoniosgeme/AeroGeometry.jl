module AeroGeometry
using LinearAlgebra
using Statistics

import Base.show

include("./Tools.jl")
include("./Units.jl")

include("./Airfoils.jl")
export Airfoil, coordinates, area,
       centroid, repanel!, write_file, local_camber, local_thickness,
       max_camber, max_thickness, leading_edge_index, trailing_edge_thickness,
       trailing_edge_angle, deflect_control_surface!,deflect_control_surface , repanel, blend_airfoils, 
       quarter_chord, list_airfoil_names

include("./Wings.jl")
export Wing, WingXSec, ControlSurface, translate!, deflect_control_surface!

include("./Fuselages.jl")
export Fuselage, FuselageXSec, coordinates, add_xsec!

include("./Airplanes.jl")
export Airplane

include("./PlotsViz.jl")

include("./Library.jl")
export cessna152

function viz end
function viz! end
export viz, viz!

end