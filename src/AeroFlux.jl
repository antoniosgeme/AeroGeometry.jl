module AeroFlux
include("./airfoil.jl")
using .airfoil

export Airfoil,plotme,get_upper_coordinates,get_lower_coordinates,get_area,
    get_centroid,repanel!,write_file,get_local_camber,get_local_thickness,
    get_max_camber,get_max_thickness, get_LE_index,get_TE_thickness
end # module AeroFlux
