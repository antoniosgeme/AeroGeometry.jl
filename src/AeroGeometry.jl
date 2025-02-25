module AeroGeometry
using LinearAlgebra
using Statistics

include("./Tools.jl")
include("./Units.jl")
export ft2m, lb2kg, psi, in2m, mi2m, mph2mps, oz2kg, gal2l, atm2pa, hp2w, kts2mps, yd2m, lbft2nm, ftlbf2j, rankine2kelvin

include("./Airfoils.jl")
export Airfoil, coordinates, area,
       centroid, repanel!, write_file, local_camber, local_thickness,
       max_camber, max_thickness, leading_edge_index, trailing_edge_thickness,
       trailing_edge_angle, deflect_control_surface!,deflect_control_surface , repanel, blend_airfoils,
       surface_coordinates, quarter_chord

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