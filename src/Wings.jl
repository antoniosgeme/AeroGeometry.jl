struct WingXSec
    airfoil::Airfoil
    xyz_le::Vector{Float64}
    chord::Float64
    twist::Float64
end 

struct Wing
    name::String
    xsecs::Vector{WingXSec}
    symmetric::Bool
end 

get_area(xsec::WingXSec) = get_area(xsec.airfoil) * xsec.chord^2
