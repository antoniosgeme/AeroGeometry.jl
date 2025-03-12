## Fuselage definitions
#==========================================================================================#

mutable struct FuselageXSec
    xyz_c::Vector{Float64}
    xyz_normal::Vector{Float64}
    width::Float64
    height::Float64
    shape::Float64

    # Custom outer constructor with keyword arguments
    function FuselageXSec(; xyz_c = [0.0, 0.0, 0.0],
                           xyz_normal = [1.0, 0.0, 0.0],
                           radius::Union{Nothing, Number} = nothing,
                           width::Union{Nothing, Number} = nothing,
                           height::Union{Nothing, Number} = nothing,
                           shape::Number = 2.0)
        if radius !== nothing
            if width !== nothing || height !== nothing
                throw(ArgumentError("Cannot specify both `radius` and (`width`, `height`) parameters - must be one or the other."))
            end
            width = 2 * radius
            height = 2 * radius
        elseif width === nothing || height === nothing
            throw(ArgumentError("Must specify either `radius` or both (`width`, `height`) parameters."))
        end

        return new(
            xyz_c,
            normalize(xyz_normal),
            width,
            height,
            shape
        )
    end
end

function show(io::IO, xsec::FuselageXSec)
    println(io, "Fuselage Cross-Section:")
    println(io, "  Center: ", xsec.xyz_c)
    println(io, "  Normal: ", xsec.xyz_normal)
    println(io, "  Width: ", xsec.width)
    println(io, "  Height: ", xsec.height)
    println(io, "  Shape: ", xsec.shape)
end


mutable struct Fuselage
    name::String
    xsecs::Vector{FuselageXSec}

    function Fuselage(;
        name::String = "Untitled",
        xsecs::Vector{FuselageXSec} = FuselageXSec[]
    )
        new(name, xsecs)
    end
end


function show(io::IO, fuselage::Fuselage)
    println(io, "Fuselage: ", fuselage.name)
    println(io, "  Number of cross-sections: ", length(fuselage.xsecs))
end


function translate!(fuselage::Fuselage, xyz::Vector{Float64})
    for xsec in fuselage.xsecs
        xsec.xyz_c .= xsec.xyz_c .+ xyz
    end 
end



function area(xsec::FuselageXSec)
    width, height, shape = xsec.width, xsec.height, xsec.shape
    return width * height / (shape^-1.8717618013591173 + 1)
end


function compute_frame(xsec::FuselageXSec)
    xyz_normal = normalize(xsec.xyz_normal)

    xg_local = xyz_normal
    zg_local = normalize([0.0, 0.0, 1.0] - dot([0.0, 0.0, 1.0], xg_local) * xg_local)
    yg_local = cross(zg_local, xg_local)

    return xg_local, yg_local, zg_local
end


function coordinates(xsec::FuselageXSec;θ::T=0:0.01:2π) where T

    st = sin.(mod.(θ, 2 * π))
    ct = cos.(mod.(θ, 2 * π))

    y = (xsec.width / 2) .* abs.(ct).^(2 / xsec.shape) .* sign.(ct)
    z = (xsec.height / 2) .* abs.(st).^(2 / xsec.shape) .* sign.(st)

    xg_local, yg_local, zg_local = compute_frame(xsec)

    x = xsec.xyz_c[1] .+ y .* yg_local[1] .+ z .* zg_local[1]
    y = xsec.xyz_c[2] .+ y .* yg_local[2] .+ z .* zg_local[2]
    z = xsec.xyz_c[3] .+ y .* yg_local[3] .+ z .* zg_local[3]

    return x, y, z
end


