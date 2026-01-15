module AeroGeometry
using LinearAlgebra
using Statistics
using Dierckx
using Printf

import Base.show

abstract type AeroComponent end

include("./Tools.jl")
include("./Units.jl")

include("./Airfoils.jl")
export Airfoil,
    coordinates,
    area,
    centroid,
    repanel!,
    write_file,
    camber,
    thickness,
    max_camber,
    max_thickness,
    leading_edge_index,
    trailing_edge_thickness,
    trailing_edge_angle,
    deflect_control_surface!,
    deflect_control_surface,
    repanel,
    blend_airfoils,
    quarter_chord,
    list_airfoil_names,
    normals,
    tangents,
    surface_coordinates,
    normalize!,
    AeroComponent

include("./ControlSurfaces.jl")
export ControlSurface

include("./Wings.jl")
export Wing, 
    WingSection, 
    translate!, 
    rotate!, 
    deflect_control_surface!, 
    span, 
    volume, 
    aerodynamic_center, 
    mean_aerodynamic_chord, 
    mean_geometric_chord

include("./Fuselages.jl")
export Fuselage, 
    FuselageSection, 
    add_xsec!, 
    ArbitraryShape,
    Hyperellipse

include("./Airplanes.jl")
export Airplane

include("Mesh.jl")
export mesh

include("./PlotsViz.jl")

include("./Library.jl")
export cessna152

function viz end
function viz! end
export viz, viz!

end
