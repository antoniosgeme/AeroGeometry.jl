module AeroGeometry
using LinearAlgebra

include("./Tools.jl")


include("./Airfoils.jl")
export Airfoil, coordinates, area,
       centroid, repanel!, write_file, local_camber, local_thickness,
       max_camber, max_thickness, leading_edge_index, trailing_edge_thickness,
       trailing_edge_angle, deflect_control_surface!,deflect_control_surface , repanel, blend_airfoils,
       surface_coordinates, quarter_chord

include("./Wings.jl")
export Wing, WingXSec, translate! 

include("./Fuselages.jl")
export Fuselage, FuselageXSec, coordinates, add_xsec!

include("./Airplanes.jl")
export Airplane


using Requires
function __init__()
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
        @info "Plots.jl detected, AeroGeometry types can now be plotted using GLMakie..."
        include("../ext/MakieVisualizations.jl")
        include("../ext/PlotsVisualizations.jl")
    end 

    @require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" begin
        @info "GLMakie.jl detected, AeroGeometry types can now be plotted using GLMakie..."
        include("../ext/MakieVisualizations.jl")
        export plotme
    end
end

end