


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


function translate!(fuselage::Fuselage, xyz::Vector{Float64})
    for xsec in fuselage.xsecs
        translate!(xsec, xyz) 
    end 
    return nothing
end

function add_xsec!(fuselage::Fuselage, xsec::FuselageXSec)
    fuselage.xsecs = vcat(fuselage.xsecs, [xsec])
    return nothing
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


function translate!(xsec::FuselageXSec, xyz::Vector{Float64})
    xsec.xyz_c = xsec.xyz_c .+ xyz,
    return nothing
end


@recipe function plot(fuselage::Fuselage; apply_limits::Bool=true)
    xlabel --> "x"
    ylabel --> "y"
    zlabel --> "z"
    legend --> :none
    markersize --> 1
    aspect_ratio --> 1
    size --> (1200, 600)
    lw --> 5

    # Deep copy fuselage to avoid modifying the original object
    fuselage_copy = deepcopy(fuselage)

    # Generate cross-section coordinates
    all_x, all_y, all_z = Float64[], Float64[], Float64[]
    for xsec in fuselage_copy.xsecs
        x, y, z = coordinates(xsec)
        append!(all_x, x)
        append!(all_y, y)
        append!(all_z, z)
    end

    # Determine limits with padding
    all_coords = vcat(all_x, all_y, all_z)
    data_min = minimum(all_coords)
    data_max = maximum(all_coords)
    padding = 0.1 * (data_max - data_min)  # Add 10% padding
    limit_min = data_min - padding
    limit_max = data_max + padding

    if apply_limits
        xlims --> (limit_min, limit_max)
        ylims --> (limit_min, limit_max)
        zlims --> (limit_min, limit_max)
    end 

    # Plot the cross-sections
    for xsec in fuselage_copy.xsecs
        x, y, z = coordinates(xsec)
        @series begin
            color := :red
            (x, y, z)
        end
    end

    # # Create a surface representation by lofting cross-sections
    if length(fuselage_copy.xsecs) > 1
        num_points = length(coordinates(fuselage_copy.xsecs[1])[1])
        num_xsecs = length(fuselage_copy.xsecs)

        # Preallocate arrays for surface coordinates
        x_surface = zeros(num_points, num_xsecs)
        y_surface = zeros(num_points, num_xsecs)
        z_surface = zeros(num_points, num_xsecs)

        for i = 1:num_xsecs
            x_surface[:, i], y_surface[:, i], z_surface[:, i] = coordinates(fuselage_copy.xsecs[i])
        end

        @series begin
            seriestype := :surface
            color := :blue          # Set color to blue
            alpha := 0.7           # Set opacity
            (x_surface, y_surface, z_surface)
        end
    end
end
