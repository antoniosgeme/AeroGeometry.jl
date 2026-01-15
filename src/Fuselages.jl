#==========================================================================================### Fuselage definitions

abstract type SectionShape end

struct Hyperellipse <: SectionShape
    width::Float64
    height::Float64
    exponent::Float64
end

struct ArbitraryShape <: SectionShape
    y::Vector{Float64}
    z::Vector{Float64}
end

mutable struct FuselageSection <: AeroComponent
    center::Vector{<:Real}
    normal::Vector{<:Real}
    shape::SectionShape
end

function FuselageSection(;
        center::Vector{<:Real} = [0.0, 0.0, 0.0],
        normal::Vector{<:Real} = [1.0, 0.0, 0.0],
        shape::Union{SectionShape, Matrix{Float64}} = Hyperellipse(2.0, 2.0, 2.0),
    )
    if shape isa Matrix{Float64}
        shape = ArbitraryShape(shape[:, 1], shape[:, 2])
    end
    return FuselageSection(center, normalize(normal), shape)
end  



function show(io::IO, section::FuselageSection)
    println(io, "Fuselage Cross-Section:")
    println(io, "  Center: ", section.center)
    println(io, "  Normal: ", section.normal)
    println(io, "  Shape: ", typeof(section.shape))
end


mutable struct Fuselage <: AeroComponent
    name::String
    sections::Vector{FuselageSection}
    position::Vector{Float64}
    orientation::Matrix{Float64}

    function Fuselage(;
        name::String = "Untitled",
        sections::Vector{FuselageSection} = FuselageSection[],
        position::Vector{Float64} = [0.0, 0.0, 0.0],
        orientation::Matrix{Float64} = Matrix{Float64}(I, 3, 3),
    )
        new(name, sections, position, orientation)
    end
end


function show(io::IO, fuselage::Fuselage)
    println(io, "Fuselage: ", fuselage.name)
    println(io, "  Number of cross-sections: ", length(fuselage.sections))
    println(io, "  Position: ", fuselage.position)
    println(io, "  Orientation:\n", fuselage.orientation)
end


function translate!(fuselage::Fuselage, xyz::Vector{<:Real})
    fuselage.position .+= xyz
    return fuselage
end


function rotate!(fuselage::Fuselage; axis::Vector{<:Real} = [1, 0, 0], angle::Real = 0.0, pivot::Vector{<:Real} = [0, 0, 0])
    pivot_global = fuselage.orientation * pivot + fuselage.position
    axis_normalized = normalize(axis)
    K = [
        0.0 -axis_normalized[3] axis_normalized[2];
        axis_normalized[3] 0.0 -axis_normalized[1];
        -axis_normalized[2] axis_normalized[1] 0.0
    ]
    R = I + sind(angle) * K + (1 - cosd(angle)) * K^2
    fuselage.position = R * (fuselage.position - pivot_global) + pivot_global
    fuselage.orientation = R * fuselage.orientation
    return fuselage
end



function area(section::FuselageSection)
    return area(section.shape)
end

function area(shape::Hyperellipse)
    width, height, shape_param = shape.width, shape.height, shape.exponent
    return width * height / (shape_param^-1.8717618013591173 + 1)
end

function area(shape::ArbitraryShape)
    x = shape.y
    y = shape.z
    return 0.5 * sum(x .* circshift(y, -1) .- y .* circshift(x, -1))

end


function compute_frame(section::FuselageSection)
    xyz_normal = normalize(section.normal)

    xg_local = xyz_normal
    zg_local = normalize([0.0, 0.0, 1.0] - dot([0.0, 0.0, 1.0], xg_local) * xg_local)
    yg_local = cross(zg_local, xg_local)

    return xg_local, yg_local, zg_local
end


function coordinates(section::FuselageSection; θ::T = 0:0.01:2π) where {T}
    y,z = coordinates(section.shape; θ = θ)

    xg_local, yg_local, zg_local = compute_frame(section)

    x = section.center[1] .+ y .* yg_local[1] .+ z .* zg_local[1]
    y = section.center[2] .+ y .* yg_local[2] .+ z .* zg_local[2]
    z = section.center[3] .+ y .* yg_local[3] .+ z .* zg_local[3]
    return x, y, z
end


function coordinates(section::Hyperellipse; θ::T = 0:0.01:2π) where {T}

    st = sin.(mod.(θ, 2 * π))
    ct = cos.(mod.(θ, 2 * π))

    y = (section.width / 2) .* abs.(ct) .^ (2 / section.exponent) .* sign.(ct)
    z = (section.height / 2) .* abs.(st) .^ (2 / section.exponent) .* sign.(st)

    return y, z
end


function coordinates(section::ArbitraryShape; θ::T = 0:0.01:2π) where {T}
    y_orig = section.y
    z_orig = section.z

    # Compute parametric angle for each point in the original shape
    θ_orig = atan.(z_orig, y_orig)
    
    # Ensure angles are in [0, 2π) and sorted
    θ_orig = mod.(θ_orig, 2π)
    sort_idx = sortperm(θ_orig)
    θ_orig = θ_orig[sort_idx]
    y_orig = y_orig[sort_idx]
    z_orig = z_orig[sort_idx]
    
    # Close the curve by appending first point at end with angle 2π
    θ_orig = vcat(θ_orig, θ_orig[1] + 2π)
    y_orig = vcat(y_orig, y_orig[1])
    z_orig = vcat(z_orig, z_orig[1])
    
    # Use Dierckx spline interpolation for smooth interpolation
    θ_target = mod.(collect(θ), 2π)
    
    spl_y = Spline1D(θ_orig, y_orig, k=3, bc="nearest")
    spl_z = Spline1D(θ_orig, z_orig, k=3, bc="nearest")
    
    y = evaluate(spl_y, θ_target)
    z = evaluate(spl_z, θ_target)

    return y, z
end

function coordinates(fuselage::Fuselage; θ::T = 0:0.01:2π) where {T}
    nθ = length(θ)
    M = length(fuselage.sections)
    X = Array{Float64}(undef, nθ, M)
    Y = Array{Float64}(undef, nθ, M)
    Z = Array{Float64}(undef, nθ, M)

    for (j, section) in enumerate(fuselage.sections)
        x_local, y_local, z_local = coordinates(section; θ = θ)
        for i in 1:nθ
            p_global = fuselage.orientation * [x_local[i]; y_local[i]; z_local[i]] + fuselage.position
            X[i, j] = p_global[1]
            Y[i, j] = p_global[2]
            Z[i, j] = p_global[3]
        end
    end

    return X, Y, Z
end



