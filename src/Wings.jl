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


    function transform_coords(coords, chord, twist_angle, le_loc)
        # Scale coordinates by the chord length
        scaled_coords = coords .* chord

        # Shift coordinates to twist about the quarter chord
        quarter_chord = 0.25 * chord
        shifted_coords = scaled_coords .- [quarter_chord 0]

        # Convert twist angle to radians
        twist_angle = deg2rad(twist_angle)

        # Rotation matrix for twisting around the spanwise (y) axis
        R = [cos(twist_angle) 0 sin(twist_angle);
                0                1 0;
            -sin(twist_angle) 0 cos(twist_angle)]

        # Apply the rotation
        rotated_coords = R * [shifted_coords[:, 1] zeros(size(shifted_coords[:, 1])) shifted_coords[:, 2]]'

        # Shift back to the quarter chord and add the leading edge location
        final_coords = rotated_coords .+ [quarter_chord + le_loc[1], le_loc[2], le_loc[3]]

        return final_coords'
    end

    # Repanel all airfoils to have the same number of points
    for xsec in wing.xsecs
        repanel!(xsec.airfoil,100)
    end

    # Generate surface data (with twist and chord applied)
    x_surface = hcat([transform_coords(xsec.airfoil.coordinates, xsec.chord, xsec.twist, xsec.le_loc)[:, 1] for xsec in wing.xsecs]...)
    y_surface = hcat([transform_coords(xsec.airfoil.coordinates, xsec.chord, xsec.twist, xsec.le_loc)[:, 2] for xsec in wing.xsecs]...)
    z_surface = hcat([transform_coords(xsec.airfoil.coordinates, xsec.chord, xsec.twist, xsec.le_loc)[:, 3] for xsec in wing.xsecs]...)

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



get_area(xsec::WingXSec) = get_area(xsec.airfoil) * xsec.chord^2


