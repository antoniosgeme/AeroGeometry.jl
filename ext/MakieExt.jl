module MakieExt

using AeroGeometry
import AeroGeometry: viz, viz! 
using Makie

@recipe(Viz, object) do scene
    Attributes(
        showsections=true,
    )
end

Makie.args_preferred_axis(::Type{<:Viz}, obj::Airplane) = LScene
Makie.args_preferred_axis(::Type{<:Viz}, obj::Wing) = LScene
Makie.args_preferred_axis(::Type{<:Viz}, obj::Fuselage) = LScene
Makie.args_preferred_axis(::Type{<:Viz}, obj::Airfoil) = Axis


"""
    Makie.plot!(plot::Viz{<:Tuple{Airplane}})

Plots an entire airplane by calling the appropriate plotting functions for its components.
"""
function Makie.plot!(plot::Viz{<:Tuple{Airplane}})
    airplane = plot[:object][]  # Extract the airplane object

    scene = parent(plot)
    scene.backgroundcolor[] = to_color(:black)  
    scene.clear[] = true  

    wing_color = RGBAf(0.1, 0.1, 0.8, 0.4)  
    fuselage_color = RGBAf(0.9, 0.5, 0.1, 0.4)  

    for fuselage in airplane.fuselages
        plot_fuselage!(plot, fuselage, fuselage_color)
    end

    for wing in airplane.wings
        plot_wing!(plot, wing, wing_color)
    end

    return plot
end

"""
    Makie.plot!(plot::Viz{<:Tuple{Wing}})

Plots only a wing.
"""
function Makie.plot!(plot::Viz{<:Tuple{Wing}})
    scene = parent(plot)
    scene.backgroundcolor[] = to_color(:black)  
    scene.clear[] = true  
    wing = plot[:object][]  
    plot_wing!(plot, wing, RGBAf(0.1, 0.1, 0.8, 0.4))
    return plot
end

"""
    Makie.plot!(plot::Viz{<:Tuple{Wing}})

Plots only a wing.
"""
function Makie.plot!(plot::Viz{<:Tuple{Vector{Wing}}})
    scene = parent(plot)
    scene.backgroundcolor[] = to_color(:black)  
    scene.clear[] = true  
    wing = plot[:object][]  
    plot_wing!(plot, wing, RGBAf(0.1, 0.1, 0.8, 0.4))
    return plot
end

"""
    Makie.plot!(plot::Viz{<:Tuple{Fuselage}})

Plots only a fuselage.
"""
function Makie.plot!(plot::Viz{<:Tuple{Fuselage}})
    scene = parent(plot)
    scene.backgroundcolor[] = to_color(:black)  
    scene.clear[] = true  
    fuselage = plot[:object][]  
    plot_fuselage!(plot, fuselage, RGBAf(0.9, 0.5, 0.1, 0.4))
    return plot
end

"""
    Makie.plot!(plot::Viz{<:Tuple{Airfoil}})

Plots only an airfoil.
"""
function Makie.plot!(plot::Viz{<:Tuple{Airfoil}})
    scene = parent(plot)
    airfoil = plot[:object][]  

    plot_airfoil!(plot, airfoil)
    return plot
end



"""
    plot_fuselage!(plot, fuselage, fuselage_color)

Plots a single fuselage, drawing its cross-sections as lines and surfaces.
"""
function plot_fuselage!(plot, fuselage, fuselage_color)
    if plot[:showsections][]
        for xsec in fuselage.xsecs
            x, y, z = coordinates(xsec)
            push!(x, first(x))
            push!(y, first(y))
            push!(z, first(z))
            lines!(plot, x, y, z, color=:red, linewidth=3)
        end
    end 

    if length(fuselage.xsecs) > 1
        num_points = length(coordinates(fuselage.xsecs[1])[1])
        num_xsecs  = length(fuselage.xsecs)

        x_surface = zeros(num_points+1, num_xsecs)
        y_surface = zeros(num_points+1, num_xsecs)
        z_surface = zeros(num_points+1, num_xsecs)

        for i in 1:num_xsecs
            xs, ys, zs = coordinates(fuselage.xsecs[i])
            x_surface[1:end-1, i] .= xs
            y_surface[1:end-1, i] .= ys
            z_surface[1:end-1, i] .= zs
            x_surface[end, i] = xs[1]
            y_surface[end, i] = ys[1]
            z_surface[end, i] = zs[1]
        end

        color_matrix = fill(fuselage_color, size(x_surface))
        surface!(plot, x_surface, y_surface, z_surface;
                 transparency=true, alpha=0.9, color=color_matrix)
    end
end

"""
    plot_wing!(plot, wing, wing_color)

Plots a single wing, drawing airfoil cross-sections and the wing surface.
"""
function plot_wing!(plot, wing, wing_color)
    x_surface, y_surface, z_surface = coordinates(wing)
    if plot[:showsections][]
        for i in 1:length(wing.xsecs)
            x_coords = x_surface[:, i]
            y_coords = y_surface[:, i]
            z_coords = z_surface[:, i]
            lines!(plot, x_coords, y_coords, z_coords, color=:red, linewidth=3)
        end
    end 

    color_matrix = fill(wing_color, size(x_surface))
    surface!(plot, x_surface, y_surface, z_surface;
             transparency=true, alpha=0.9, color=color_matrix)

    if wing.symmetric
        if plot[:showsections][]
            for i in 1:length(wing.xsecs)
                x_coords = x_surface[:, i]
                y_coords = y_surface[:, i]
                z_coords = z_surface[:, i]
                lines!(plot, x_coords, -y_coords, z_coords, color=:red, linewidth=3)
            end
        end 
        surface!(plot, x_surface, -y_surface, z_surface;
                 invert_normals=true, transparency=true, alpha=0.9, color=color_matrix)
    end
end

"""
    plot_airfoil!(plot, airfoil)

Plots an airfoil using its coordinates.
"""
function plot_airfoil!(plot, airfoil)
    xy = coordinates(airfoil)
    lines!(plot, xy[:,1], xy[:,2], color=:blue, linewidth=2)
end

end # module


