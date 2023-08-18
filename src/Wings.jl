mutable struct WingXSec
    airfoil::Airfoil
    xyz_le::Vector{Float64}
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

get_area(xsec::WingXSec) = get_area(xsec.airfoil) * xsec.chord^2
