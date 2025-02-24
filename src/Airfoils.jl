## Fuselage definitions
#==========================================================================================#

using Dierckx
using Printf

mutable struct Airfoil
    name::String
    coordinates::Matrix{Float64}      # The current (possibly deflected) coordinates
    original_coordinates::Matrix{Float64}  # The original, undeflected coordinates
end

Airfoil(name::String, coordinates::Matrix{Float64}) = Airfoil(name, coordinates, coordinates)

"""
    Airfoil(name::String)

Constructs an Airfoil object from the provided name. The function follows these steps to generate the airfoil:
1. If the name is a valid file path, the coordinates are loaded from the file.
2. If the name corresponds to a NACA 4-digit airfoil, its coordinates are generated directly.
3. If neither of the above succeeds, an attempt is made to find the airfoil in the UIUC database.

The Airfoil object contains the name of the airfoil and a 2D array of `x` and `y` coordinates, starting from the trailing edge (TE), moving along the upper surface to the leading edge (LE), and then along the lower surface back to the TE.

# Examples
```julia-repl
julia> hub = Airfoil("naca0012")  # Generates a NACA 0012 airfoil
julia> tip = Airfoil("b707d")     # Finds and loads an airfoil from the UIUC database
julia> custom = Airfoil("my_airfoil.dat")  # Loads airfoil coordinates from a file
```
"""
function Airfoil(name::String)

    # Check if the file exists and load coordinates from it
    if isfile(name) || isfile(name * ".dat")
        coordinates = from_filepath(name)
        if !isnothing(coordinates)
            return Airfoil(name, coordinates)
        end
    end
    # Handle NACA airfoils
    if occursin("naca", lowercase(name)) && length(filter(isdigit, name)) == 4
        coordinates = naca_coords(name)
        if !isnothing(coordinates)
            return Airfoil(name, coordinates)
        end
    end

    # If not found, check in the UIUC database
    coordinates = UIUC_coords(name)
    if !isnothing(coordinates)
        return Airfoil(name, coordinates)
    end

    println("Unable to generate airfoil. Returning empty Airfoil object...")
    return Airfoil("None",zeros(1,1))
end

"""
Internal function used to populate airfoil coordinates from a file path.
"""
function from_filepath(name::String)
    # Ensure the filename ends with ".dat"
    datfile = strip(name)
    if !endswith(datfile, ".dat")
        datfile *= ".dat"
    end

    # Read the file and extract coordinates
    try
        io = readlines(datfile)
        x = [parse(Float64, split(line)[1]) for line in io[2:end]]
        y = [parse(Float64, split(line)[2]) for line in io[2:end]]
        return hcat(x, y)
    catch e
        println("Error reading file \"$datfile\": $e")
        return nothing
    end
end

"""
Internal function used to populate airfoil coordinates from UIUC database.
"""
function UIUC_coords(name)
    # Ensure the filename ends with ".dat"
    datfile = strip(lowercase(name))
    if !endswith(datfile, ".dat")
        datfile *= ".dat"
    end
    
    # Construct the directory path
    dir = joinpath(dirname(@__FILE__), "airfoil_database")
    
    if datfile in readdir(dir)
        filepath = joinpath(dir, datfile)
        try
            io = readlines(filepath)
            x = [parse(Float64, split(xx)[1]) for xx in io[2:end]]
            y = [parse(Float64, split(yy)[2]) for yy in io[2:end]]
            return hcat(x, y)
        catch
            println("Error reading the airfoil file.")
            return nothing
        end
    else
        println("Airfoil file not found in UIUC database.")
        return nothing
    end
end


"""
    naca_coords(name::String, points_per_side::Int64 = 50)

Generates coordinates for 4 digit NACA airfoils with 2*points_per_side-1 points 
    using cosine spacing
"""
function naca_coords(name::String, points_per_side::Int64 = 100)

    naca_num = strip(name)[5:end]

    if length(naca_num) != 4
        Error("Oops! Can only populate from 4 digit naca airfoils")
    end

    m = parse(Float64,string(naca_num[1]))/100.
    p = parse(Float64,string(naca_num[2]))/10.
    t = parse(Float64,string(naca_num[3:end]))/100.

    x_t = cos_space(0,1,points_per_side) 

    x_t1 = x_t[x_t .<= p]
    x_t2 = x_t[x_t .> p]


    y_t(x) = 5t*(0.2969 * √x - 0.1260x - 0.3516x^2 + 0.2843x^3 - 0.1015x^4) # 0.1015/0.1036 for blunt/sharp TE

    p == 0  ? p = 0.5 : nothing

    y_c1(x) =(m / p^2) * (2p * (x) - x^2)
    y_c2(x) = m / (1 - p)^2 * ((1 - 2p) + 2p * x - x^2)

    y_c = vcat(y_c1.(x_t1) , y_c1.(x_t2))

    dyc1_dx(x) = 2m / p^2 * (p-x)
    dyc2_dx(x) = 2m / (1-p)^2 * (p-x)

    dyc_dc = vcat(dyc1_dx.(x_t1),dyc2_dx.(x_t2))

    θ = atan.(dyc_dc)

    x_U = similar(x_t)
    x_L = similar(x_t)
    y_U = similar(x_t)
    y_L = similar(x_t)

    @. x_U = x_t - y_t(x_t) * sin(θ)
    @. x_L = x_t + y_t(x_t) * sin(θ)
    @. y_U = y_c + y_t(x_t) * cos(θ)
    @. y_L = y_c - y_t(x_t) * cos(θ)

    # Flip upper surface so it's back to front
    x_U = x_U[end:-1:1, :]
    y_U = y_U[end:-1:1, :]

    # Trim 1 point from lower surface so there's no overlap
    x_L = x_L[2:end]
    y_L = y_L[2:end]

    x = vcat(x_U, x_L)
    y = vcat(y_U, y_L)


    return hcat(x,y)
end



"""
    leading_edge_index(airfoil::Airfoil)

Retrieves index of the leading edge in the coordinate array of an Airfoil object
"""
leading_edge_index(airfoil::Airfoil) = argmin(airfoil.coordinates[:,1])

"""
    max_camber(airfoil::Airfoil;x_over_c=0:0.01:1)

Retrieves the maximum camber an Airfoil object
"""
max_camber(airfoil::Airfoil;x_over_c=0:0.01:1) = maximum(local_camber(airfoil,x_over_c=x_over_c))

"""
    max_thickness(airfoil::Airfoil;x_over_c=0:0.01:1)

Retrieves the maximum thickness an Airfoil object
"""
max_thickness(airfoil::Airfoil;x_over_c=0:0.01:1) = maximum(local_thickness(airfoil,x_over_c=x_over_c))


"""
    coordinates(airfoil::Airfoil, surface::Symbol=:all)

Retrieves the coordinates of an Airfoil object based on the specified surface.

Arguments:

- `airfoil`: The Airfoil object.
- `surface`: A Symbol indicating which surface to retrieve. Options are:
  - `:upper`: Retrieves the upper surface of the airfoil.
  - `:lower`: Retrieves the lower surface of the airfoil.
  - `:all` (default): Retrieves all coordinates.

Returns:

A matrix containing the coordinates of the specified surface.
"""
function coordinates(airfoil::Airfoil, surface::Symbol=:all)::Matrix{Float64}
    split_index = argmin(airfoil.coordinates[:, 1])

    if surface == :upper
        return airfoil.coordinates[split_index:-1:1, :]
    elseif surface == :lower
        return airfoil.coordinates[split_index:end, :]
    elseif surface == :all
        return airfoil.coordinates
    else
        throw(ArgumentError("Invalid surface option. Choose from :upper, :lower, or :all."))
    end
end

"""
    area(airfoil::Airfoil)

Retrieves the area an Airfoil object
"""
area(airfoil::Airfoil) = 0.5 * sum(airfoil.coordinates[:,1] .* circshift(airfoil.coordinates[:,2],-1)
                                 .- airfoil.coordinates[:,2] .* circshift(airfoil.coordinates[:,1],-1))

 
"""
    surface_coordinates(airfoil::Airfoil)

Retrieves the coordinate along the surface of an Airfoil object, starting at the trailing edge 
of the upper surface and wrapping around the leading edge to end at the trailing edge 
of the lower surface
"""                                 
function surface_coordinates(airfoil::Airfoil)
    s = zeros(size(airfoil.coordinates,1))
    dx = diff(airfoil.coordinates[:,1])
    dy = diff(airfoil.coordinates[:,2])
    ds = hypot.(dx,dy)
    s[2:end] = cumsum(ds)
    return s
end 


"""
    local_camber(airfoil::Airfoil;x_over_c=0:0.01:1)

Retrieves the camber distribution of an Airfoil object
""" 
function local_camber(airfoil::Airfoil;x_over_c=0:0.01:1)
    upper = coordinates(airfoil,:upper)
    lower = coordinates(airfoil,:lower)
    interpolator_lower = Spline1D(lower[:,1],lower[:,2],bc="nearest")
    interpolator_upper = Spline1D(upper[:,1],upper[:,2],bc="nearest")
    interp_upper = evaluate(interpolator_upper,x_over_c)
    interp_lower = evaluate(interpolator_lower,x_over_c)
    return ( interp_upper + interp_lower )/2
end 

"""
    local_thickness(airfoil::Airfoil;x_over_c=0:0.01:1)

Retrieves the thickness distribution of an Airfoil object
""" 
function local_thickness(airfoil::Airfoil;x_over_c=0:0.01:1)
    upper = coordinates(airfoil,:upper)
    lower = coordinates(airfoil,:lower)
    interpolator_lower = Spline1D(lower[:,1],lower[:,2],bc="nearest")
    interpolator_upper = Spline1D(upper[:,1],upper[:,2],bc="nearest")
    interp_upper = evaluate(interpolator_upper,x_over_c)
    interp_lower = evaluate(interpolator_lower,x_over_c)
    return ( interp_upper - interp_lower ) 
end 


"""
    trailing_edge_thickness(airfoil::Airfoil)

Retrieves the trailing edge thickness of an Airfoil object
""" 
function trailing_edge_thickness(airfoil::Airfoil)
    x_gap = airfoil.coordinates[1,1] - airfoil.coordinates[end,1]
    y_gap = airfoil.coordinates[1,2] - airfoil.coordinates[end,2]
    return hypot(x_gap,y_gap)
end 

"""
    trailing_edge_angle(airfoil::Airfoil)

Retrieves the trailing edge angle of an Airfoil object
""" 
function trailing_edge_angle(airfoil)
    upper_vec = airfoil.coordinates[1,:] - airfoil.coordinates[2,:]
    lower_vec = airfoil.coordinates[end,:] - airfoil.coordinates[end-1,:]
    return atand(
        upper_vec[1] * lower_vec[2] - upper_vec[2] * lower_vec[1],
        upper_vec[1] * lower_vec[1] + upper_vec[2] * upper_vec[2]
        )
end 

"""
    centroid(airfoil::Airfoil)

Retrieves the geometric centroid of an Airfoil object
""" 
function centroid(airfoil::Airfoil)
    x =  airfoil.coordinates[:,1]
    y  = airfoil.coordinates[:,2]

    xn = circshift(airfoil.coordinates[:,1],-1)
    yn = circshift(airfoil.coordinates[:,2],-1)

    a = x .* yn .- y .* xn

    A = 0.5 * sum(a)

    xc = 1 / (6 * A) * sum(a .* (x .+ xn))
    yc = 1 / (6 * A) * sum(a .* (y .+ yn))

    return xc , yc
end


"""
quarter_chord(airfoil::Airfoil)

Retrieves the geometric quarter chord point of an Airfoil object
""" 
function quarter_chord(airfoil::Airfoil)
    # Find the leading edge (minimum x-coordinate)
    le_index = leading_edge_index(airfoil)
    leading_edge = airfoil.coordinates[le_index, :]

    # Find the trailing edge (maximum x-coordinate)
    te_index = argmax(airfoil.coordinates[:,1])
    trailing_edge = airfoil.coordinates[te_index, :]

    # Compute quarter chord point
    quarter_chord_x = leading_edge[1] + 0.25 * (trailing_edge[1] - leading_edge[1])
    quarter_chord_y = leading_edge[2] + 0.25 * (trailing_edge[2] - leading_edge[2])

    return [quarter_chord_x, quarter_chord_y]
end




"""
    repanel!(airfoil::Airfoil,points_per_side)  

Repanels an Airfoil object in place according to points\\_per\\_side. The total
number of points will be (2*points\\_per\\_side - 1)
""" 
function repanel!(airfoil::Airfoil, points_per_side; hinge=nothing)
    LE_index = leading_edge_index(airfoil)
    s = surface_coordinates(airfoil)
    x = airfoil.coordinates[:,1]
    y = airfoil.coordinates[:,2]

    x_interpolator = Spline1D(s, x, bc="nearest")
    y_interpolator = Spline1D(s, y, bc="nearest")

    if hinge isa Number && 0 < hinge < 1
        # Convert hinge fraction to arc-length coordinate
        s_hinge = minimum(s) + hinge * (maximum(s) - minimum(s))

        # Find the closest index to the hinge
        hinge_index = searchsortedfirst(s, s_hinge)

        # Define consistent numbers of points
        points_before_hinge = div(points_per_side, 2)
        points_after_hinge = points_per_side - points_before_hinge

        # Generate new panel distributions
        s_upper_new = vcat(
            cos_space(minimum(s[1:LE_index]), s_hinge, points_before_hinge),
            cos_space(s_hinge, maximum(s[1:LE_index]), points_after_hinge)
        )

        s_lower_new = vcat(
            cos_space(minimum(s[LE_index:end]), s_hinge, points_before_hinge),
            cos_space(s_hinge, maximum(s[LE_index:end]), points_after_hinge)
        )

    else
        # Regular cosine spacing if no hinge is specified
        s_upper_new = cos_space(minimum(s[1:LE_index]), maximum(s[1:LE_index]), points_per_side)
        s_lower_new = cos_space(minimum(s[LE_index:end]), maximum(s[LE_index:end]), points_per_side)
    end

    # Combine upper and lower surfaces
    s_new = vcat(s_upper_new[1:end-1], s_lower_new)
    x_new = evaluate(x_interpolator, s_new)
    y_new = evaluate(y_interpolator, s_new)

    # Update airfoil coordinates
    airfoil.coordinates = hcat(x_new, y_new)

    return airfoil
end


"""
    repanel(airfoil::Airfoil,points_per_side)  

Create a new Airfoil object, repanelled according to points\\_per\\_side. The total
number of points will be (2*points\\_per\\_side - 1)
""" 
function repanel(airfoil::Airfoil,points_per_side)
    airfoil_new = deepcopy(airfoil)
    repanel!(airfoil_new,points_per_side)
    return airfoil_new
end

"""
    write_file(airfoil::Airfoil)  

Creates a .dat file of an Airfoil object for use in other software
""" 
function write_file(airfoil::Airfoil)
    open(airfoil.name*".dat","w") do f
        write(f,airfoil.name*"\n")
        for i in 1:length(airfoil.coordinates[:,1])
            var = (airfoil.coordinates[i,1],airfoil.coordinates[i,2])
            str = "$(@sprintf("     %.12f    %.12f\n",var...))"
            write(f,str)
        end
    end
end


"""
    deflect_control_surface!(airfoil::Airfoil; deflection=0, x_hinge=0.75)

Adds a control surface the trailing edge of an Airfoil object. deflection 
specifies the deflection angle in degrees, and x_hinge specifies how far down 
the chord the hinge is to be located 
""" 
function deflect_control_surface!(airfoil::Airfoil; deflection=0, x_hinge=0.75)
    # Compute hinge point
    camber = local_camber(airfoil, x_over_c=x_hinge)
    thickness = local_thickness(airfoil, x_over_c=x_hinge)
    y_hinge = camber + (thickness / 2) * (deflection <= 0 ? 1 : -1)
    hinge_point = [x_hinge, y_hinge]

    # Retrieve coordinates
    upper = airfoil.original_coordinates[argmin(airfoil.original_coordinates[:,1]):-1:1,:]
    lower = airfoil.original_coordinates[argmin(airfoil.original_coordinates[:,1]):end,:]

    # Define rotation function
    function rotate_and_translate(coords, hinge_point, angle)
        delta = coords .- hinge_point'
        rotated = rotate2D(delta, -angle) .+ hinge_point'
        return rotated
    end

    # Rotate points
    upper_rotated = rotate_and_translate(upper, hinge_point, deflection)
    lower_rotated = rotate_and_translate(lower, hinge_point, deflection)

    # Identify points near the hinge
    smoothing_radius = 0.1  # Region around hinge for smoothing
    upper_blend_mask = abs.(upper[:, 1] .- x_hinge) .< smoothing_radius
    lower_blend_mask = abs.(lower[:, 1] .- x_hinge) .< smoothing_radius

    # Smooth points near hinge using Dierckx
    function smooth_points(x, y, mask)
        spline = Spline1D(x[mask], y[mask], k=3)  # Cubic spline
        smoothed_y = spline(x[mask])
        return smoothed_y
    end

    upper[upper_blend_mask, 2] .= smooth_points(
        upper[:, 1], upper[:, 2], upper_blend_mask
    )
    upper[upper_blend_mask, 2] .= smooth_points(
        upper[:, 1], upper[:, 2], upper_blend_mask
    )

    # Update rotated points behind the hinge
    upper_behind = upper[:, 1] .>= x_hinge
    lower_behind = lower[:, 1] .>= x_hinge

    upper[upper_behind, :] .= upper_rotated[upper_behind, :]
    lower[lower_behind, :] .= lower_rotated[lower_behind, :]

    # Update airfoil coordinates
    new_x = vcat(reverse(upper[:, 1]), lower[2:end, 1])
    new_y = vcat(reverse(upper[:, 2]), lower[2:end, 2])
    airfoil.coordinates = hcat(new_x, new_y)

    return airfoil
end

"""
    deflect_control_surface(airfoil::Airfoil; deflection=0, x_hinge=0.75)

Returns an Airfoil object with a control surface at the trailing edge. deflection 
specifies the deflection angle in degrees, and x_hinge specifies how far down 
the chord the hinge is to be located 
""" 
function deflect_control_surface(airfoil::Airfoil; deflection=0, x_hinge=0.75)
    airfoil_new = deepcopy(airfoil)
    deflect_control_surface!(airfoil_new,deflection=deflection,x_hinge=x_hinge)
    return airfoil_new
end 


"""
    blend_airfoils(airfoil1::Airfoil,airfoil2::Airfoil;fraction::Number=0.5,points_per_side=100)

Superimposes two airfoils, taking fraction of the first airfoil and 1-fraction of the second. 
Repanels according to points\\_per\\_side. See repanel! for more info
""" 
function blend_airfoils(airfoil1::Airfoil,airfoil2::Airfoil;fraction::Number=0.5,points_per_side=100)
    repaneled1 = repanel(airfoil1,points_per_side)
    repaneled2 = repanel(airfoil2,points_per_side)
    coordinates = repaneled1.coordinates * fraction + repaneled2.coordinates * (1 - fraction)
    return Airfoil(airfoil1.name*"+"*airfoil2.name,coordinates)
end 
