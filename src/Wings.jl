module Wings

export WingXSec, Wing

using ..Airfoils

struct WingXSec
    airfoil::Airfoil
    xyz_le::Vector{Float64}
    chord::Float64
    twist::Float64
end 

#get_area(xsec::WingXSec) = get_area(xsec.airfoil) * xsec.chord^2


struct Wing
    name::String
    xsecs::Vector{WingXSec}
    symmetric::Bool
end 

end 