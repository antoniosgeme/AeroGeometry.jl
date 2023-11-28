using RecipesBase

mutable struct WingXSec
    airfoil::Airfoil
    le_loc::Vector{Float64}
    chord::Float64
    twist::Float64
end 

WingXSec(airfoil::Airfoil) = WingXSec(airfoil,[0,0,0],1,0)
WingXSec(airfoil::Airfoil,chord::Float64) = WingXSec(airfoil,[0,0,0],chord,0)

mutable struct Wing
    name::String
    xsecs::Vector{WingXSec}
    symmetric::Bool
end 

@recipe function plot(wing::Wing) 
    xlabel --> "x" 
    ylabel --> "y"
    zlabel --> "z"
    zlims --> (-2,2)
    markersize --> 1
    fillrange --> 0
    aspect_ratio --> 1
    legend --> :none
    for i in eachindex(wing.xsecs)
        @series begin
            (
                wing.xsecs[i].airfoil.coordinates[:,1] .+ wing.xsecs[i].le_loc[1], 
                zeros(size(wing.xsecs[i].airfoil.coordinates[:,2])) .+ wing.xsecs[i].le_loc[2],
                wing.xsecs[i].airfoil.coordinates[:,2] .+ wing.xsecs[i].le_loc[3]
            )
        end
    end 
end


get_area(xsec::WingXSec) = get_area(xsec.airfoil) * xsec.chord^2


