module Wings

export WingXSec

using ..Airfoils

struct WingXSec
    airfoil::Airfoil
    xyz_le::Vector{Float64}
    chord::Float64
    twist::Float64
end 
end 