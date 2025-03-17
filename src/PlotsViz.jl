using RecipesBase
@recipe function plot(airfoil::Airfoil; camberline=false)
    xlabel --> "x"
    ylabel --> "y"
    markersize --> 1
    aspect_ratio --> 1
    legend --> :none
    size --> (1600, 800)
    bg --> :black
    lw --> 3
    grid --> true  # Enable grid
    gridlinewidth --> 3
    gridstyle --> :dash
    gridcolor --> :white  # Set white grid lines
    
    # Adjust minor grid visibility
    minorgrid --> true
    minorgridcolor --> :white
    minorgridlinewidth --> 0.4
    minorgridstyle --> :dot

    @series begin
        (
            airfoil.coordinates[:,1], 
            airfoil.coordinates[:,2]
        )
    end
    if camberline
        @series begin
            lw := 1.5
            linestyle := :dot  # Differentiate camberline
            (
                0:0.01:1, 
                local_camber(airfoil)
            )
        end
    end 
end

@recipe function plot(wing::Wing;isometric=true)
    xlabel --> "x [m]"
    ylabel --> "y [m]"
    zlabel --> "z [m]"
    legend --> :none
    marksize --> 1
    size --> (1600, 800)
    lw --> 3
    bg --> :black

    if isometric
        airplane = Airplane(wings=[wing])
        xm, ym, zm, d = isometric_limits(airplane)
        xlims --> (xm-d,xm+d)
        ylims --> (ym-d,ym+d)
        zlims --> (zm-d,zm+d)
    end 

    wing_copy = deepcopy(wing)
    # Repanel all airfoils to have the same number of points
    for xsec in wing_copy.xsecs
        repanel!(xsec.airfoil,100)
    end

    # Generate surface data (with twist and chord applied)
    (x_surface,y_surface,z_surface) = coordinates(wing_copy)


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



@recipe function plot(fuselage::Fuselage; isometric=true)
    xlabel --> "x [m]"
    ylabel --> "y [m]"
    zlabel --> "z [m]"
    legend --> :none
    size --> (1600, 800)
    lw --> 3
    bg --> :black

    if isometric
        airplane = Airplane(fuselages=[fuselage])
        xm, ym, zm, d = isometric_limits(airplane)
        xlims --> (xm-d,xm+d)
        ylims --> (ym-d,ym+d)
        zlims --> (zm-d,zm+d)
    end 

    

    # Deep copy fuselage to avoid modifying the original object
    fuselage_copy = deepcopy(fuselage)

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





@recipe function plot(airplane::Airplane; isometric=true)
    xlabel --> "x [m]"
    ylabel --> "y [m]"
    zlabel --> "z [m]"
    legend --> :none
    size --> (1600, 800)
    markersize --> 1
    lw --> 3
    bg --> :black
    title --> airplane.name

    if isometric
        xm, ym, zm, d = isometric_limits(airplane)
        xlims --> (xm-d,xm+d)
        ylims --> (ym-d,ym+d)
        zlims --> (zm-d,zm+d)
    end 

    # Plot all fuselages
    for fuselage in airplane.fuselages
        for xsec in fuselage.xsecs
            x, y, z = coordinates(xsec)
            push!(x, first(x))  # Close the loop for cross-section
            push!(y, first(y))
            push!(z, first(z))
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
            x_surface = zeros(num_points+1, num_xsecs)
            y_surface = zeros(num_points+1, num_xsecs)
            z_surface = zeros(num_points+1, num_xsecs)

            for i = 1:num_xsecs
                x_surface[1:end-1, i], y_surface[1:end-1, i], z_surface[1:end-1, i] = coordinates(fuselage.xsecs[i])
                x_surface[end, i] = x_surface[1, i]  # Close the loop
                y_surface[end, i] = y_surface[1, i]
                z_surface[end, i] = z_surface[1, i]
            end

            @series begin
                seriestype := :surface
                color := :orange          # Set color to blue
                alpha := 0.7           # Set opacity
                (x_surface, y_surface, z_surface)
            end
        end
    end 



    # plot wings
    for wing in airplane.wings
        wing_copy = deepcopy(wing)
        # Repanel all airfoils to have the same number of points
        #for xsec in wing_copy.xsecs
            #repanel!(xsec.airfoil,100)
        #end

    # Generate surface data (with twist and chord applied)
        (x_surface,y_surface,z_surface) = coordinates(wing_copy)
        for i = axes(x_surface,2)
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



 # Initialize arrays to track all coordinates
function isometric_limits(airplane::Airplane)
    all_x, all_y, all_z = Float64[], Float64[], Float64[]
    for fuse in airplane.fuselages
        for xsec in fuse.xsecs
            x, y, z = coordinates(xsec)
            append!(all_x, x)
            append!(all_y, y)
            append!(all_z, z)
        end
    end
    for wing in airplane.wings
        x_surface, y_surface, z_surface = coordinates(wing)
        append!(all_x, x_surface...)
        append!(all_y, y_surface...)
        append!(all_z, z_surface...)
        if wing.symmetric
            y_surface_symmetric = -y_surface
            append!(all_x, x_surface...)
            append!(all_y, y_surface_symmetric...)
            append!(all_z, z_surface...)
        end
    end
    x12, y12, z12 = extrema(all_x), extrema(all_y), extrema(all_z)
    d = maximum([diff([x12...]),diff([y12...]),diff([z12...])])[1] / 2
    xm, ym, zm = mean(x12),  mean(y12),  mean(z12) 
    return xm, ym, zm, d
end 