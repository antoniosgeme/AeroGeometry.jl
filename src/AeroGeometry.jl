module AeroGeometry
using LinearAlgebra
using Statistics

import Base.show

abstract type AeroComponent end

include("./tools.jl")
include("./units.jl")

include("./airfoils.jl")
export Airfoil, coordinates, area,
       centroid, repanel!, write_file, camber, thickness,
       max_camber, max_thickness, leading_edge_index, trailing_edge_thickness,
       trailing_edge_angle, deflect_control_surface!,deflect_control_surface , repanel, blend_airfoils, 
       quarter_chord, list_airfoil_names, normals, tangents, surface_coordinates, AeroComponent

include("./wings.jl")
export Wing, WingSection, ControlSurface, translate!, deflect_control_surface!

include("./fuselages.jl")
export Fuselage, FuselageSection, coordinates, add_xsec!

include("./airplanes.jl")
export Airplane

include("./plotsviz.jl")

include("./library.jl")
export cessna152

function viz end
function viz! end
export viz, viz!

end