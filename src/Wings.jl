using RecipesBase

mutable struct WingXSec
    airfoil::Airfoil
    le_loc::Vector{Float64}
    chord::Float64
    twist::Float64
end 

WingXSec(airfoil::Airfoil) = WingXSec(airfoil,[0,0,0],1,0)
WingXSec(airfoil::Airfoil,chord::Float64) = WingXSec(airfoil,[0,0,0],chord,0)

    mutable struct Wing
        name::String
        xsecs::Vector{WingXSec}
        symmetric::Bool
    end 

@recipe function plot(wing::Wing)
    xlabel --> "x"
    ylabel --> "y"
    zlabel --> "z"
    legend --> :none
    markersize --> 1
    aspect_ratio --> 1
    size --> (1200, 600)
    lw --> 3


    # Repanel all airfoils to have the same number of points
    for xsec in wing.xsecs
        repanel!(xsec.airfoil,100)
    end

    # Generate surface data (with twist and chord applied)
    (x_surface,y_surface,z_surface) = get_global_coordinates(wing::Wing)

    # Determine limits with padding
    all_coords = vcat(x_surface, y_surface, z_surface)
    data_min = minimum(all_coords)
    data_max = maximum(all_coords)
    padding = 0.1 * (data_max - data_min)  # Add 10% padding
    limit_min = data_min - padding
    limit_max = data_max + padding

    if wing.symmetric
        limit_min = -limit_max
    end 
    xlims --> (limit_min, limit_max)
    ylims --> (limit_min, limit_max)
    zlims --> (limit_min, limit_max)

     # Plot the cross-sections
     for i = 1:length(wing.xsecs)
        x_coords = x_surface[:, i]
        y_coords = y_surface[:, i]
        z_coords = z_surface[:, i]
        @series begin
            color := :red
            (
                x_coords, 
                y_coords,
                z_coords
            )
        end

        if wing.symmetric && i>1
            x_coords = x_surface[:, i]
            y_coords = y_surface[:,1]-y_surface[:, i]
            z_coords = z_surface[:, i]
            @series begin
                color := :red
                (
                    x_coords, 
                    y_coords,
                    z_coords
                )
            end
        end 
    end

    #Plot the surface
    @series begin
        seriestype := :surface
        color := :blue          # Set color to blue
        alpha := 0.7           # Set opacity
        (x_surface, y_surface, z_surface)
    end

    if wing.symmetric
        @series begin
            seriestype := :surface
            color := :blue          # Set color to blue
            alpha := 0.7           # Set opacity
            (x_surface, y_surface[:,1].-y_surface, z_surface)
        end
    end 
end


"""
    get_global_coordinates(wing::Wing)

Computes the coordinates of each wing cross-section in the global reference frame
"""
function get_global_coordinates(wing::Wing)
    N = length(wing.xsecs[1].airfoil.coordinates[:,1])
    x_surface = zeros(N,length(wing.xsecs))
    y_surface = zeros(N,length(wing.xsecs))
    z_surface = zeros(N,length(wing.xsecs))

    for i in 1:length(wing.xsecs)
        xsec = wing.xsecs[i]
        coords = hcat(xsec.airfoil.coordinates[:,1], zeros(size(xsec.airfoil.coordinates[:,1])), xsec.airfoil.coordinates[:,2])
        chord = xsec.chord
        twist_angle = xsec.twist
        le_loc = xsec.le_loc

        # Scale coordinates by the chord length
        scaled_coords = coords .* chord

        # Move the leading edge to the `le_loc`
        translated_coords = scaled_coords .+ le_loc'

        if i < length(wing.xsecs)
            next_xsec = wing.xsecs[i + 1]
            direction = 1
        else
            next_xsec = wing.xsecs[i - 1]
            direction = -1
        end 
        # Find the axis of rotation
        qc = get_quarter_chord(xsec.airfoil) .* chord
        qc_next = get_quarter_chord(next_xsec.airfoil) .* chord
        quarter_chord_current = le_loc .+ [qc[1],0,qc[2]]
        quarter_chord_next = next_xsec.le_loc .+ [qc_next[1],0,qc_next[2]]
        axis = direction * normalize(quarter_chord_next - quarter_chord_current)

        # Convert twist angle to radians
        twist_angle_rad = deg2rad(twist_angle)

        # Rotate each point around the computed axis
        for (j, point) in enumerate(eachrow(translated_coords))
            rotated_coords = rotate_vector(point - quarter_chord_current, axis, twist_angle_rad) + quarter_chord_current
            x_surface[j,i] = rotated_coords[1]
            y_surface[j,i] = rotated_coords[2]
            z_surface[j,i] = rotated_coords[3]
        end 
    end
    return (x_surface, y_surface, z_surface)
end



get_area(xsec::WingXSec) = get_area(xsec.airfoil) * xsec.chord^2


