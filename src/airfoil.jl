module airfoil

export Airfoil,plotme,get_upper_coordinates,get_lower_coordinates,get_area,
    get_centroid,repanel!,write_file,get_local_camber,get_local_thickness,
    get_max_camber,get_max_thickness, get_LE_index,get_TE_thickness,
    get_TE_angle,add_control_surface!,repanel


include(".\\tools.jl")

using .tools

using Plots
using PCHIPInterpolation
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
    dir = replace(@__FILE__,"airfoil.jl" =>"airfoil_database\\")
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


function plotme(foil::Airfoil)

    fig = scatter(foil.coordinates[:,1],foil.coordinates[:,2],markersize=1)
    scatter!(aspect_ratio = :equal,xlabel = "x",ylabel = "y",legend=false)
    return fig
end

function plotme(arr::Array{<:Number,2})
    fig = scatter(arr[:,1],arr[:,2])
    scatter!(aspect_ratio = :equal,xlabel = "x",ylabel = "y",legend=false,markersize=1)
    return fig
end

function plotme(x::Array{<:Number,1},y::Array{<:Number,1})
    fig = scatter(x,y)
    scatter!(aspect_ratio = :equal,xlabel = "x",ylabel = "y",legend=false,markersize=1)
    return fig
end

get_upper_coordinates(foil::Airfoil) = foil.coordinates[argmin(foil.coordinates[:,1]):-1:1,:]
get_upper_coordinates(coordinates::Array{<:Number,2}) = coordinates[argmin(coordinates[:,1]):-1:1,:]

get_lower_coordinates(foil::Airfoil) = foil.coordinates[argmin(foil.coordinates[:,1]):end,:]
get_lower_coordinates(coordinates::Array{<:Number,2}) = coordinates[argmin(coordinates[:,1]):end,:]


get_area(foil::Airfoil) = 0.5 * sum(foil.coordinates[:,1] .* circshift(foil.coordinates[:,2],-1)
                                 .- foil.coordinates[:,2] .* circshift(foil.coordinates[:,1],-1))
                                

function get_local_camber(foil::Airfoil;x_over_c=0:0.01:1)
    upper = get_upper_coordinates(foil)
    lower = get_lower_coordinates(foil)
    interpolator_lower = Interpolator(lower[:,1],lower[:,2])
    interpolator_upper = Interpolator(upper[:,1],upper[:,2])
    interp_upper = interpolator_upper.(x_over_c)
    interp_lower = interpolator_lower.(x_over_c)
    return ( interp_upper + interp_lower )/2
end 

function get_local_thickness(foil::Airfoil;x_over_c=0:0.01:1)
    upper = get_upper_coordinates(foil)
    lower = get_lower_coordinates(foil)
    interpolator_lower = Interpolator(lower[:,1],lower[:,2])
    interpolator_upper = Interpolator(upper[:,1],upper[:,2])
    interp_upper = interpolator_upper.(x_over_c)
    interp_lower = interpolator_lower.(x_over_c)
    return ( interp_upper - interp_lower ) 
end 

get_LE_index(foil) = argmin(foil.coordinates[:,1])
get_max_camber(foil,;x_over_c=0:0.01:1) = maximum(get_local_camber(foil,x_over_c=x_over_c))
get_max_thickness(foil,;x_over_c=0:0.01:1) = maximum(get_local_thickness(foil,x_over_c=x_over_c))

function get_TE_thickness(foil)
    x_gap = foil.coordinates[1,1] - foil.coordinates[end,1]
    y_gap = foil.coordinates[1,2] - foil.coordinates[end,2]
    return hypot(x_gap,y_gap)
end 

function get_TE_angle(foil)
    upper_vec = foil.coordinates[1,:] - foil.coordinates[2,:]
    lower_vec = foil.coordinates[end,:] - foil.coordinates[end-1,:]
    return atand(
        upper_vec[1] * lower_vec[2] - upper_vec[2] * lower_vec[1],
        upper_vec[1] * lower_vec[1] + upper_vec[2] * upper_vec[2]
        )
end 

function get_centroid(foil::Airfoil)
    x =  foil.coordinates[:,1]
    y  = foil.coordinates[:,2]

    xn = circshift(foil.coordinates[:,1],-1)
    yn = circshift(foil.coordinates[:,2],-1)

    a = x .* yn .- y .* xn

    A = 0.5 * sum(a)

    x_c = 1 / (6 * A) * sum(a .* (x .+ xn))
    y_c = 1 / (6 * A) * sum(a .* (y .+ yn))

    return x_c , y_c
end

function repanel(coordinates::Array{<:Number,2},points_per_side)
    upper = get_upper_coordinates(coordinates)
    lower = get_lower_coordinates(coordinates)
    interp_upper = Interpolator(upper[:,1],upper[:,2])
    interp_lower = Interpolator(lower[:,1],lower[:,2])
    s_upper = cos_space(minimum(upper[:,1])+eps(),maximum(upper[:,1])-eps(),points_per_side)
    s_lower = cos_space(minimum(lower[:,1])+eps(),maximum(lower[:,1])-eps(),points_per_side)
    new_upper = interp_upper.(s_upper)
    new_lower = interp_lower.(s_lower)
    new_x = vcat(s_upper[end:-1:1],s_lower[2:end])
    new_y = vcat(new_upper[end:-1:1],new_lower[2:end])
    return hcat(new_x,new_y)
end 

function repanel!(foil::Airfoil,points_per_side)
    foil.coordinates = repanel(foil.coordinates,points_per_side)
    return nothing
end

function repanel(foil::Airfoil,points_per_side)
    coordinates= repanel(foil.coordinates,points_per_side)
    return Airfoil(foil.name,coordinates)
end

function write_file(foil::Airfoil)
    open(foil.name*".dat","w") do f
        write(f,foil.name*"\n")
        for i in 1:length(foil.coordinates[:,1])
            var = (foil.coordinates[i,1],foil.coordinates[i,2])
            str = "$(@sprintf("     %.12f    %.12f\n",var...))"
            write(f,str)
        end
    end
end



function add_control_surface!(foil::Airfoil; deflection=0, x_hinge=0.75)
    if deflection > 0
        y_hinge =  get_local_camber(foil,x_over_c=x_hinge) - get_local_thickness(foil,x_over_c=x_hinge) / 2
    else
        y_hinge =  get_local_camber(foil,x_over_c=x_hinge) + get_local_thickness(foil,x_over_c=x_hinge) / 2
    end 

    hinge_point = [x_hinge,y_hinge]
    upper = get_upper_coordinates(foil)
    lower = get_lower_coordinates(foil)

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
    foil.coordinates = hcat(new_x,new_y)
    return nothing
end 

function blend_airfoils(foil1::Airfoil,foil2::Airfoil,fraction)


end 


end
