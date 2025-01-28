import GLMakie

function plotme(airplane::Airplane)
    wing_color = GLMakie.RGBAf(0.1, 0.1, 0.8, 0.4)  
    fuselage_color = GLMakie.RGBAf(0.9, 0.5, 0.1, 0.4)  # Soft orange

    # Create a Figure and LScene for interactivity
    fig = GLMakie.Figure(size = (1200, 800))
    pl = GLMakie.PointLight(GLMakie.Point3f(0, 0, 10), GLMakie.RGBf(1.0, 1.0, 1.0))  # Neutral point light
    al = GLMakie.AmbientLight(GLMakie.RGBf(1, 1, 1))
    lscene = GLMakie.LScene(fig[1, 1], show_axis=false, scenekw = (lights = [pl, al], backgroundcolor=:black, clear=true))

    # Plot fuselages
    for fuselage in airplane.fuselages
        
        for xsec in fuselage.xsecs
            x, y, z = coordinates(xsec)
            push!(x, first(x))  # Close the loop for cross-section
            push!(y, first(y))
            push!(z, first(z))
            GLMakie.lines!(lscene, x, y, z, color=:red, linewidth=3)
        end

        # Loft surface for fuselage if it has multiple cross-sections
        if length(fuselage.xsecs) > 1
            num_points = length(coordinates(fuselage.xsecs[1])[1])
            num_xsecs = length(fuselage.xsecs)

            x_surface = zeros(num_points+1, num_xsecs)
            y_surface = zeros(num_points+1, num_xsecs)
            z_surface = zeros(num_points+1, num_xsecs)

            for i = 1:num_xsecs
                x_surface[1:end-1, i], y_surface[1:end-1, i], z_surface[1:end-1, i] = coordinates(fuselage.xsecs[i])
                x_surface[end, i] = x_surface[1, i]  # Close the loop
                y_surface[end, i] = y_surface[1, i]
                z_surface[end, i] = z_surface[1, i]
            end

            color_matrix = [ fuselage_color for i in 1:size(x_surface, 1), j in 1:size(x_surface, 2)]
            GLMakie.surface!(lscene, x_surface, y_surface, z_surface, transparency=true, alpha=0.9,color=color_matrix)
        end
    end

    # Plot wings
    for wing in airplane.wings
        # Ensure consistent airfoil resolution
        wing_copy = deepcopy(wing)
        for xsec in wing_copy.xsecs
            repanel!(xsec.airfoil, 100)
        end

        x_surface, y_surface, z_surface = coordinates(wing_copy)

        # Plot each section of the wing as lines
        for i = 1:length(wing_copy.xsecs)
            x_coords = x_surface[:, i]
            y_coords = y_surface[:, i]
            z_coords = z_surface[:, i]
            GLMakie.lines!(lscene, x_coords, y_coords, z_coords, color=:red, linewidth=3)
        end

        # Plot the wing surface
        
        color_matrix = [ wing_color for i in 1:size(x_surface, 1), j in 1:size(x_surface, 2)]
        GLMakie.surface!(lscene, x_surface, y_surface, z_surface,color=color_matrix,transparency=true, alpha=0.9)

        # Plot symmetric wing surface if applicable
        if wing.symmetric
            for i = 1:length(wing_copy.xsecs)
                x_coords = x_surface[:, i]
                y_coords = y_surface[:, i]
                z_coords = z_surface[:, i]
                GLMakie.lines!(lscene, x_coords, -y_coords, z_coords, color=:red, linewidth=3)
            end
            GLMakie.surface!(lscene, x_surface, -y_surface, z_surface,invert_normals=true,transparency=true, alpha=0.9,color=color_matrix)
        end
    end

    return fig
end

