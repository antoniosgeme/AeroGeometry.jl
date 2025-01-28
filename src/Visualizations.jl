using RecipesBase

@recipe function plot(airfoil::Airfoil) 
    xlabel --> "x" 
    ylabel --> "y"
    markersize --> 1
    aspect_ratio --> 1
    legend --> :none
    return airfoil.coordinates[:,1],airfoil.coordinates[:,2]
end

@recipe function plot(wing::Wing ; apply_limits::Bool=true)
    xlabel --> "x"
    ylabel --> "y"
    zlabel --> "z"
    legend --> :none
    markersize --> 1
    aspect_ratio --> 1
    size --> (1200, 600)
    lw --> 3

    wing_copy = deepcopy(wing)
    # Repanel all airfoils to have the same number of points
    for xsec in wing_copy.xsecs
        repanel!(xsec.airfoil,100)
    end

    # Generate surface data (with twist and chord applied)
    (x_surface,y_surface,z_surface) = coordinates(wing_copy)

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
    if apply_limits
        xlims --> (limit_min, limit_max)
        ylims --> (limit_min, limit_max)
        zlims --> (limit_min, limit_max)
    end 

    # Plot the cross-sections
    for i = 1:length(wing_copy.xsecs)
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

        if wing_copy.symmetric && i>1
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



@recipe function plot(airplane::Airplane)
    xlabel --> "x"
    ylabel --> "y"
    zlabel --> "z"
    legend --> :none
    markersize --> 1
    aspect_ratio --> 1
    size --> (1600, 800)
    lw --> 3

    # Initialize arrays to track all coordinates
    all_x, all_y, all_z = Float64[], Float64[], Float64[]

    # Process fuselages
    for fuse in airplane.fuselages
        for xsec in fuse.xsecs
            x, y, z = coordinates(xsec)
            append!(all_x, x)
            append!(all_y, y)
            append!(all_z, z)
        end
    end

    # Process wings
    for wing in airplane.wings
        x_surface, y_surface, z_surface = coordinates(wing)
        append!(all_x, x_surface...)
        append!(all_y, y_surface...)
        append!(all_z, z_surface...)

        # Handle symmetry for wings
        if wing.symmetric
            y_surface_symmetric = -y_surface
            append!(all_x, x_surface...)
            append!(all_y, y_surface_symmetric...)
            append!(all_z, z_surface...)
        end
    end


    data_min = minimum(vcat(all_x, all_y, all_z))
    data_max = maximum(vcat(all_x, all_y, all_z))

    padding = 0.1 * (data_max - data_min)  # Add 10% padding

    

    xlims --> (data_min-padding, data_max+padding)
    ylims --> (data_min-padding, data_max+padding)
    zlims --> (data_min-padding, data_max+padding)

    # Plot all fuselages
    for fuselage in airplane.fuselages
        for xsec in fuselage.xsecs
            x, y, z = coordinates(xsec)
            @series begin
                color := :red
                (x, y, z)
            end
        end

        # # Create a surface representation by lofting cross-sections
        if length(fuselage.xsecs) > 1
            num_points = length(coordinates(fuselage.xsecs[1])[1])
            num_xsecs = length(fuselage.xsecs)

            # Preallocate arrays for surface coordinates
            x_surface = zeros(num_points, num_xsecs)
            y_surface = zeros(num_points, num_xsecs)
            z_surface = zeros(num_points, num_xsecs)

            for i = 1:num_xsecs
                x_surface[:, i], y_surface[:, i], z_surface[:, i] = coordinates(fuselage.xsecs[i])
            end

            @series begin
                seriestype := :surface
                color := :blue          # Set color to blue
                alpha := 0.7           # Set opacity
                (x_surface, y_surface, z_surface)
            end
        end
    end 



    # plot wings
    for wing in airplane.wings
        wing_copy = deepcopy(wing)
        # Repanel all airfoils to have the same number of points
        for xsec in wing_copy.xsecs
            repanel!(xsec.airfoil,100)
        end

    # Generate surface data (with twist and chord applied)
    (   x_surface,y_surface,z_surface) = coordinates(wing_copy)
        for i = 1:length(wing_copy.xsecs)
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

            if wing_copy.symmetric && i>1
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
end