using Revise
using AeroGeometry
import GLMakie

#import Plots as plt
#plt.plotly()

ft(feet, inches) = 0.3048 * feet + 0.0254 * inches

# Define the wings
wing_xsecs = [
    WingXSec(
        airfoil=Airfoil("naca2412"),
        le_loc=[0, 0, 0],
        chord=ft(5,4),
    ),
    WingXSec(
        airfoil=Airfoil("naca2412"),
        le_loc=[0, ft(7,0), ft(7,0) * sind(1)],
        chord=ft(5,4),
    ),
    WingXSec(
        airfoil=Airfoil("naca0012"),
        le_loc=[ft(4, 3/4) - ft(3, 8 + 1/2), ft(33, 4)/2, ft(33, 4)/2 * sind(1)],
        chord=ft(3, 8 + 1/2),
        twist=0
    )
]
wing = Wing(name="Main Wing", xsecs=wing_xsecs, symmetric=true)

hs_xsecs = [
    WingXSec(
        airfoil=Airfoil("naca0012"),
        le_loc=[0, 0, 0],
        chord=ft(3,8),
        twist=-2
    ),
    WingXSec(
        airfoil=Airfoil("naca0012"),
        le_loc=[ft(1,0), ft(10,0) / 2, 0],
        chord=ft(2, 4 + 3 / 8),
        twist=-2
    )
]
horizontal_stabilizer = Wing(name="Horizontal Stabilizer", xsecs=hs_xsecs, symmetric=true)
translate!(horizontal_stabilizer, [4.0648, 0, -0.6096])

vs_xsecs = [
    WingXSec(
        airfoil=Airfoil("naca0012"),
        le_loc=[ft(-5,0), 0, 0],
        chord=ft(8, 8),
        twist=0
    ),
    WingXSec(
        airfoil=Airfoil("naca0012"),
        le_loc=[ft(0,0), 0, ft(1,0)],
        chord=ft(3, 8),
        twist=0
    ),
    WingXSec(
        airfoil=Airfoil("naca0012"),
        le_loc=[ft(0, 8), 0, ft(5,0)],
        chord=ft(2, 8), 
        twist=0
    )
]
vertical_stabilizer = Wing(name="Vertical Stabilizer", xsecs=vs_xsecs)
translate!(vertical_stabilizer, [ft(16, 11) - ft(3, 8), 0, ft(-2,0)])

# Define the fuselage
xc =[0,0,ft(3,0),ft(5,0),ft(10,4),ft(21,11)]
zc = [ft(-1,0),ft(-1,0),ft(-0.85,0),ft(0,0),ft(0.3,0),ft(0.8,0)]
radii = [0, ft(1.5,0), ft(1.7,0), ft(2.7,0), ft(2.3,0), ft(0.3,0)]
shapes = [2, 3, 7, 7, 7, 3]

fuse_xsecs = [FuselageXSec(radius=radii[i],xyz_c=[xc[i], 0, -0.3048],shape=shapes[i]) for i in eachindex(xc)]
fuselage = Fuselage(name="Main Body", xsecs=fuse_xsecs)
translate!(fuselage, [ft(-5,0), 0, ft(-2,0)])
# Combine into an airplane
airplane = Airplane(
    name="Cessna 152",
    wings=[wing, horizontal_stabilizer, vertical_stabilizer],
    fuselages=[fuselage]
)



function plot_airplane_with_lscene(airplane::Airplane)
    # Create a Figure and LScene for interactivity
    fig = GLMakie.Figure(size = (1200, 800))
    pl = GLMakie.PointLight(GLMakie.Point3f(0, 0, 10), GLMakie.RGBf(1.0, 1.0, 1.0))  # Neutral point light
    al = GLMakie.AmbientLight(GLMakie.RGBf(1, 1, 1))
    lscene = GLMakie.LScene(fig[1, 1], show_axis=false, scenekw = (lights = [pl, al], backgroundcolor=:white, clear=true))

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

            #GLMakie.surface!(lscene, x_surface, y_surface, z_surface, color=:blue, transparency=true, alpha=0.7,shading=GLMakie.NoShading)
            
            # Wireframe for fuselage loft
            #for i = 1:num_xsecs
                #GLMakie.lines!(lscene, x_surface[:, i], y_surface[:, i], z_surface[:, i], color=:blue, linewidth=1)
            #end
            for j = 1:10:num_points+1
                GLMakie.lines!(lscene, x_surface[j, :], y_surface[j, :], z_surface[j, :], color=:gray, linewidth=1)
            end
        end
    end

   # Plot wings as wireframes
   for wing in airplane.wings
    # Ensure consistent airfoil resolution
    wing_copy = deepcopy(wing)
    for xsec in wing_copy.xsecs
        repanel!(xsec.airfoil, 100)
    end

    x_surface, y_surface, z_surface = coordinates(wing_copy)

    # Plot wireframe for the wing
    for i = 1:size(x_surface, 2)
        GLMakie.lines!(lscene, x_surface[:, i], y_surface[:, i], z_surface[:, i], color=:red, linewidth=1)
    end
    for j = 1:size(x_surface, 1)
        GLMakie.lines!(lscene, x_surface[j, :], y_surface[j, :], z_surface[j, :], color=:red, linewidth=1)
    end

    # Symmetric wing
    if wing.symmetric
        for i = 1:10:size(x_surface, 2)
            GLMakie.lines!(lscene, x_surface[:, i], -y_surface[:, i], z_surface[:, i], color=:red, linewidth=1)
        end
        for j = 1:10:size(x_surface, 1)
            GLMakie.lines!(lscene, x_surface[j, :], -y_surface[j, :], z_surface[j, :], color=:red, linewidth=1)
        end
    end
end

return fig
end

# Plot the airplane
plot_airplane_with_lscene(airplane)

