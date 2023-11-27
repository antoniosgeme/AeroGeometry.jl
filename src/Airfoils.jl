include("./Tools.jl")
using .Tools

using RecipesBase
using Dierckx
using Printf

mutable struct Airfoil
    name::String
    coordinates::Array{Float64,2}
end

"""
    Airfoil(name::String ; filepath::String = "")

Construccts an Airfoil object from filepath provided. If filepath is empty
an attempt to find the airfoil in the UIUC database is made and if that fails
an attempt to generate it as a naca 4 digit airfoil is made

Airfoil object contains the name of Airfoil and a 2D Array of x and y
coordinates starting from TE and ending at TE.

# Example
```julia-repl
julia> hub = Airfoil("naca0012")
julia> tip = Airfoil("b707d")
```
"""
function Airfoil(name::String ; filepath::String = "")
    if ! isempty(filepath)
        coordinates = from_filepath(filepath)
        if ! isnothing(coordinates)
            return Airfoil(name,coordinates)
        else
            println("File not found...")
        end 
    else
        if occursin("naca",lowercase(name)) && length(filter(isdigit,name))==4
            coordinates = naca_coords(name)
            return Airfoil(name,coordinates)   
        end 

        coordinates = UIUC_coords(name)
        if ! isnothing(coordinates)
            return Airfoil(name,coordinates)
        end

        println("Unable to generate airfoil...")
        return nothing
    end
end

"""
Internal function used to populate airfoil coordinates from fileapath.
"""
function from_filepath(name)

    if  ! occursin(".dat",name)
        datfile = strip(name)*".dat"
    else
        datfile = strip(name)
    end

    try
        f = open(datfile,"r")
    catch
        return nothing
    end 

    io = readlines(f)

    close(f)

    x = [parse(Float64,split(xx)[1]) for xx in io[2:end]]
    y = [parse(Float64,split(yy)[2]) for yy in io[2:end]]

    return hcat(x,y)
end

"""
Internal function used to populate airfoil coordinates from UIUC database.
"""
function UIUC_coords(name)

    datfile = strip(lowercase(name))*".dat"
    dir = replace(@__FILE__,"Airfoils.jl" =>"airfoil_database\\")
    if datfile in readdir(dir)
        f = open(dir*datfile,"r")
        io = readlines(f)
        close(f)

        x = [parse(Float64,split(xx)[1]) for xx in io[2:end]]
        y = [parse(Float64,split(yy)[2]) for yy in io[2:end]]

        return hcat(x,y)

    else
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

@recipe function plot(airfoil::Airfoil) 
    xlabel --> "x" 
    ylabel --> "y"
    markersize --> 1
    fillrange --> 0
    aspect_ratio --> 1
    legend --> :none
    return airfoil.coordinates[:,1],airfoil.coordinates[:,2]
end

get_upper_coordinates(airfoil::Airfoil) = airfoil.coordinates[argmin(airfoil.coordinates[:,1]):-1:1,:]
get_upper_coordinates(coordinates::Array{<:Number,2}) = coordinates[argmin(coordinates[:,1]):-1:1,:]

get_lower_coordinates(airfoil::Airfoil) = airfoil.coordinates[argmin(airfoil.coordinates[:,1]):end,:]
get_lower_coordinates(coordinates::Array{<:Number,2}) = coordinates[argmin(coordinates[:,1]):end,:]


get_area(airfoil::Airfoil) = 0.5 * sum(airfoil.coordinates[:,1] .* circshift(airfoil.coordinates[:,2],-1)
                                 .- airfoil.coordinates[:,2] .* circshift(airfoil.coordinates[:,1],-1))

function get_surface_coordinates(airfoil::Airfoil)
    s = zeros(size(airfoil.coordinates,1))
    dx = diff(airfoil.coordinates[:,1])
    dy = diff(airfoil.coordinates[:,2])
    ds = hypot.(dx,dy)
    s[2:end] = cumsum(ds)
    return s
end 

function get_surface_coordinates(coordinates::Array{<:Number,2})
    s = zeros(size(coordinates,1))
    dx = diff(coordinates[:,1])
    dy = diff(coordinates[:,2])
    ds = hypot.(dx,dy)
    s[2:end] = cumsum(ds)
    return s
end 

function get_local_camber(airfoil::Airfoil;x_over_c=0:0.01:1)
    upper = get_upper_coordinates(airfoil)
    lower = get_lower_coordinates(airfoil)
    interpolator_lower = Spline1D(lower[:,1],lower[:,2],bc="nearest")
    interpolator_upper = Spline1D(upper[:,1],upper[:,2],bc="nearest")
    interp_upper = evaluate(interpolator_upper,x_over_c)
    interp_lower = evaluate(interpolator_lower,x_over_c)
    return ( interp_upper + interp_lower )/2
end 

function get_local_thickness(airfoil::Airfoil;x_over_c=0:0.01:1)
    upper = get_upper_coordinates(airfoil)
    lower = get_lower_coordinates(airfoil)
    interpolator_lower = Spline1D(lower[:,1],lower[:,2],bc="nearest")
    interpolator_upper = Spline1D(upper[:,1],upper[:,2],bc="nearest")
    interp_upper = evaluate(interpolator_upper,x_over_c)
    interp_lower = evaluate(interpolator_lower,x_over_c)
    return ( interp_upper - interp_lower ) 
end 

get_LE_index(airfoil::Airfoil) = argmin(airfoil.coordinates[:,1])
get_LE_index(coordinates::Array{<:Number,2}) = argmin(coordinates[:,1])

get_max_camber(airfoil::Airfoil,;x_over_c=0:0.01:1) = maximum(get_local_camber(airfoil,x_over_c=x_over_c))
get_max_thickness(airfoil::Airfoil,;x_over_c=0:0.01:1) = maximum(get_local_thickness(airfoil,x_over_c=x_over_c))

function get_TE_thickness(airfoil)
    x_gap = airfoil.coordinates[1,1] - airfoil.coordinates[end,1]
    y_gap = airfoil.coordinates[1,2] - airfoil.coordinates[end,2]
    return hypot(x_gap,y_gap)
end 

function get_TE_angle(airfoil)
    upper_vec = airfoil.coordinates[1,:] - airfoil.coordinates[2,:]
    lower_vec = airfoil.coordinates[end,:] - airfoil.coordinates[end-1,:]
    return atand(
        upper_vec[1] * lower_vec[2] - upper_vec[2] * lower_vec[1],
        upper_vec[1] * lower_vec[1] + upper_vec[2] * upper_vec[2]
        )
end 

function get_centroid(airfoil::Airfoil)
    x =  airfoil.coordinates[:,1]
    y  = airfoil.coordinates[:,2]

    xn = circshift(airfoil.coordinates[:,1],-1)
    yn = circshift(airfoil.coordinates[:,2],-1)

    a = x .* yn .- y .* xn

    A = 0.5 * sum(a)

    x_c = 1 / (6 * A) * sum(a .* (x .+ xn))
    y_c = 1 / (6 * A) * sum(a .* (y .+ yn))

    return x_c , y_c
end

function repanel(coordinates::Array{<:Number,2},points_per_side)
    LE_index = get_LE_index(coordinates)
    s = get_surface_coordinates(coordinates)
    s_upper = s[1:LE_index]
    s_lower = s[LE_index:end]
    x = coordinates[:,1]
    y = coordinates[:,2]
    x_interpolator = Spline1D(s,x,bc="nearest")
    y_interpolator = Spline1D(s,y,bc="nearest")
    s_upper_new = cos_space(minimum(s_upper),maximum(s_upper),points_per_side)
    s_lower_new = cos_space(minimum(s_lower),maximum(s_lower),points_per_side)
    s_new = vcat(s_upper_new[1:end-1],s_lower_new)
    x_new = evaluate(x_interpolator,s_new)
    y_new = evaluate(y_interpolator,s_new)
    return hcat(x_new,y_new)
end 

function repanel!(airfoil::Airfoil,points_per_side)
    airfoil.coordinates = repanel(airfoil.coordinates,points_per_side)
    return nothing
end

function repanel(airfoil::Airfoil,points_per_side)
    coordinates= repanel(airfoil.coordinates,points_per_side)
    return Airfoil(airfoil.name,coordinates)
end

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



function add_control_surface!(airfoil::Airfoil; deflection=0, x_hinge=0.75)
    if deflection > 0
        y_hinge =  get_local_camber(airfoil,x_over_c=x_hinge) - get_local_thickness(airfoil,x_over_c=x_hinge) / 2
    else
        y_hinge =  get_local_camber(airfoil,x_over_c=x_hinge) + get_local_thickness(airfoil,x_over_c=x_hinge) / 2
    end 

    hinge_point = [x_hinge,y_hinge]
    upper = get_upper_coordinates(airfoil)
    lower = get_lower_coordinates(airfoil)

    upper_rotated = zeros(size(upper))
    lower_rotated = zeros(size(lower))
    for i in eachindex(upper[:,1])
        upper_rotated[i,:] = rotate2D(upper[i,:] .- hinge_point, -deflection) .+ hinge_point  
    end 
    for i in eachindex(lower[:,1])
        lower_rotated[i,:] = rotate2D(lower[i,:] .- hinge_point, -deflection) .+ hinge_point 
    end 

    upper_behind = ( ( upper_rotated[:,1] .- x_hinge ) * cosd(deflection/2)  - 
                     ( upper_rotated[:,2] .- y_hinge ) * sind(deflection/2) ) .>= 0
                     
    lower_behind = ( ( lower_rotated[:,1] .- x_hinge ) * cosd(deflection/2) - 
                     ( lower_rotated[:,2] .- y_hinge ) * sind(deflection/2) ) .>= 0

    upper[upper_behind,:] = upper_rotated[upper_behind,:]
    lower[lower_behind,:] = lower_rotated[lower_behind,:]

    new_x = vcat(upper[end:-1:1,1],lower[2:end,1])
    new_y = vcat(upper[end:-1:1,2],lower[2:end,2])
    airfoil.coordinates = hcat(new_x,new_y)
    return nothing
end 

function blend_airfoils(airfoil1::Airfoil,airfoil2::Airfoil;fraction::Number=0.5,points_per_side=100)
    repaneled1 = repanel(airfoil1,points_per_side)
    repaneled2 = repanel(airfoil2,points_per_side)
    coordinates = repaneled1.coordinates * fraction + repaneled2.coordinates * (1 - fraction)
    return Airfoil(airfoil1.name*"+"*airfoil2.name,coordinates)
end 
