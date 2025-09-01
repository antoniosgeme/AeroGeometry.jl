## Fuselage definitions
#==========================================================================================#

mutable struct FuselageSection <: AeroComponent
    center::Vector{Float64}
    normal::Vector{Float64}
    width::Float64
    height::Float64
    shape::Float64

    # Custom outer constructor with keyword arguments
    function FuselageSection(; center = [0.0, 0.0, 0.0],
                           normal = [1.0, 0.0, 0.0],
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
            center,
            normalize(normal),
            width,
            height,
            shape
        )
    end
end

function show(io::IO, section::FuselageSection)
    println(io, "Fuselage Cross-Section:")
    println(io, "  Center: ", section.center)
    println(io, "  Normal: ", section.normal)
    println(io, "  Width: ", section.width)
    println(io, "  Height: ", section.height)
    println(io, "  Shape: ", section.shape)
end


mutable struct Fuselage <: AeroComponent
    name::String
    sections::Vector{FuselageSection}

    function Fuselage(;
        name::String = "Untitled",
        sections::Vector{FuselageSection} = FuselageSection[]
    )
        new(name, sections)
    end
end


function show(io::IO, fuselage::Fuselage)
    println(io, "Fuselage: ", fuselage.name)
    println(io, "  Number of cross-sections: ", length(fuselage.sections))
end


function translate!(fuselage::Fuselage, xyz::Vector{Float64})
    for xsec in fuselage.sections
        xsec.center .= xsec.center .+ xyz
    end 
end



function area(section::FuselageSection)
    width, height, shape = section.width, section.height, section.shape
    return width * height / (shape^-1.8717618013591173 + 1)
end


function compute_frame(section::FuselageSection)
    xyz_normal = normalize(section.normal)

    xg_local = xyz_normal
    zg_local = normalize([0.0, 0.0, 1.0] - dot([0.0, 0.0, 1.0], xg_local) * xg_local)
    yg_local = cross(zg_local, xg_local)

    return xg_local, yg_local, zg_local
end


function coordinates(section::FuselageSection;θ::T=0:0.01:2π) where T

    st = sin.(mod.(θ, 2 * π))
    ct = cos.(mod.(θ, 2 * π))

    y = (section.width / 2) .* abs.(ct).^(2 / section.shape) .* sign.(ct)
    z = (section.height / 2) .* abs.(st).^(2 / section.shape) .* sign.(st)

    xg_local, yg_local, zg_local = compute_frame(section)

    x = section.center[1] .+ y .* yg_local[1] .+ z .* zg_local[1]
    y = section.center[2] .+ y .* yg_local[2] .+ z .* zg_local[2]
    z = section.center[3] .+ y .* yg_local[3] .+ z .* zg_local[3]

    return x, y, z
end


