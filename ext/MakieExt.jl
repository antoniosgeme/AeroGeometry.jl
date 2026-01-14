module MakieExt

using AeroGeometry
import AeroGeometry: viz, viz!, isometric_limits
using Makie


Makie.convert_arguments(P::PointBased,airfoil::Airfoil) = convert_arguments(P,airfoil.x, airfoil.y)
Makie.convert_arguments(S::Type{<:Surface},wing::Wing) = convert_arguments(S,coordinates(wing)...)



function viz(airfoil::Airfoil; show_camber=false, show_thickness=false)

    fig, ax, plt = lines(airfoil,color=:blue,linewidth=2,label="Airfoil")
    ax.aspect = DataAspect()
    ax.xlabel = "X"
    ax.ylabel = "Y"
    ax.title = airfoil.name
    xlims!(ax, -0.1, 1.1)
    ylims!(ax, -0.6, 0.6)

    xy = coordinates(airfoil)
    fill_color=(:lightblue, 0.9)
    poly!(ax, Point2f.(xy[:, 1], xy[:, 2]), color=fill_color)
    

    if show_camber
        camb = camber(airfoil)
        lines!(ax, camb[:, 1], camb[:, 2], color=:red, linewidth=2,label="Camber line")
    end

    if show_thickness
        thick = thickness(airfoil)
        lines!(ax, thick[:, 1], thick[:, 2], color=:green, linewidth=2,label="Thickness line")
    end


    axislegend(ax)
    display(fig)
    return fig, ax, plt
end


function plot_wing!(ax, wing::Wing; show_sections=true)
    X, Y, Z = coordinates(wing)
    if show_sections
        for i = axes(X, 2)
            x_coords = X[:, i]
            y_coords = Y[:, i]
            z_coords = Z[:, i]
            lines!(ax, x_coords, y_coords, z_coords, color = :red, linewidth = 3)
        end
    end

    color_matrix = fill(:blue, size(X))
    plt = surface!(
        ax,
        X,
        Y,
        Z;
        transparency = true,
        alpha = 0.5,
        color = color_matrix,
    )
    return plt
end


function plot_fuselage!(ax, fuselage::Fuselage; show_sections=true)
    X, Y, Z = coordinates(fuselage)
    X = vcat(X, X[1, :]')
    Y = vcat(Y, Y[1, :]')
    Z = vcat(Z, Z[1, :]')
    if show_sections
        for i = axes(X, 2)
            x_coords = X[:, i]
            y_coords = Y[:, i]
            z_coords = Z[:, i]
            lines!(ax, x_coords, y_coords, z_coords, color = :red, linewidth = 3)
        end
    end

    color_matrix = fill(RGBAf(0.9, 0.5, 0.1, 0.4), size(X))
    plt = surface!(
        ax,
        X,
        Y,
        Z;
        transparency = true,
        alpha = 0.9,
        color = color_matrix,
    )
    return plt
end

function viz(wing::Wing; show_sections=true)
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect=:data,clip=false)
    plt = plot_wing!(ax, wing; show_sections=show_sections)
    airplane = Airplane(wings = [wing])
    xm, ym, zm, d = isometric_limits(airplane)
    d = 1.1d
    xlims!(ax, xm - d, xm + d)
    ylims!(ax, ym - d, ym + d)
    zlims!(ax, zm - d, zm + d)
    display(fig)
    return fig, ax, plt
end 



function viz(fuselage::Fuselage; show_sections=true)
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect=:data,clip=false)
    plt = plot_fuselage!(ax, fuselage; show_sections=show_sections)

    airplane = Airplane(fuselages = [fuselage])
    xm, ym, zm, d = isometric_limits(airplane)
    d = 1.1d
    xlims!(ax, xm - d, xm + d)
    ylims!(ax, ym - d, ym + d)
    zlims!(ax, zm - d, zm + d)
    display(fig)
    return fig, ax, plt

end 


function viz(airplane::Airplane; show_sections=true)

    fig = Figure()
    ax = Axis3(fig[1, 1], aspect=:data,clip=false)
    plt = nothing

    for wing in airplane.wings
        plot_wing!(ax, wing; show_sections=show_sections)
    end

    for fuselage in airplane.fuselages
        plt = plot_fuselage!(ax, fuselage; show_sections=show_sections)
    end

    xm, ym, zm, d = isometric_limits(airplane)
    d = 1.1d
    xlims!(ax, xm - d, xm + d)
    ylims!(ax, ym - d, ym + d)
    zlims!(ax, zm - d, zm + d)
    display(fig)
    return fig, ax, plt
end 

end # module
