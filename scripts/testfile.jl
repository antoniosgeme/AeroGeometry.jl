using AeroGeometry
import Plots as plt
plt.plotly()

import AeroGeometry: ft2m
# Define the wings
wing_sections = [
    WingSection(airfoil = Airfoil("naca2412"), le_loc = [0, 0, 0], chord = ft2m(5, 4)),
    WingSection(
        airfoil = Airfoil("naca2412"),
        le_loc = [0, ft2m(7, 0), ft2m(7, 0) * sind(1)],
        chord = ft2m(5, 4),
    ),
    WingSection(
        airfoil = Airfoil("naca0012"),
        le_loc = [ft2m(4, 3/4) - ft2m(3, 8 + 1/2), ft2m(32, 4)/2, ft2m(33, 4)/2 * sind(1)],
        chord = ft2m(3, 8 + 1/2),
        twist = 0,
    ),
]

aileron = ControlSurface(name = "Aileron", xsec_id = [2, 3], symmetric = false)
wing = Wing(
    name = "Main Wing",
    sections = wing_sections,
    symmetric = true,
    control_surfaces = [aileron],
)

hs_sections = [
    WingSection(
        airfoil = Airfoil("naca0012"),
        le_loc = [0, 0, 0],
        chord = ft2m(3, 8),
        twist = -2,
    ),
    WingSection(
        airfoil = Airfoil("naca0012"),
        le_loc = [ft2m(1, 0), ft2m(10, 0) / 2, 0],
        chord = ft2m(2, 4 + 3 / 8),
        twist = -2,
    ),
]

elevator = ControlSurface(name = "Elevator", xsec_id = [1, 2], symmetric = true)
horizontal_stabilizer = Wing(
    name = "Horizontal Stabilizer",
    sections = hs_sections,
    symmetric = true,
    control_surfaces = [elevator],
)
AeroGeometry.translate!(horizontal_stabilizer, [4.0648, 0, -0.6096])

vs_sections = [
    WingSection(
        airfoil = Airfoil("naca0012"),
        le_loc = [ft2m(-5, 0), 0, 0],
        chord = ft2m(8, 8),
        twist = 0,
    ),
    WingSection(
        airfoil = Airfoil("naca0012"),
        le_loc = [ft2m(0, 0), 0, ft2m(1, 0)],
        chord = ft2m(3, 8),
        twist = 0,
    ),
    WingSection(
        airfoil = Airfoil("naca0012"),
        le_loc = [ft2m(0, 8), 0, ft2m(5, 0)],
        chord = ft2m(2, 8),
        twist = 0,
    ),
]

rudder = ControlSurface(name = "Rudder", xsec_id = [2, 3], symmetric = true)
vertical_stabilizer = Wing(
    name = "Vertical Stabilizer",
    sections = vs_sections,
    symmetric = false,
    control_surfaces = [rudder],
)
AeroGeometry.translate!(vertical_stabilizer, [ft2m(16, 11) - ft2m(3, 8), 0, ft2m(-2, 0)])

#define the fuselage
xc = [0, 0, ft2m(3, 0), ft2m(5, 0), ft2m(10, 4), ft2m(12, 4), ft2m(21, 11)]
zc = [
    ft2m(-1, 0),
    ft2m(-1, 0),
    ft2m(-0.85, 0),
    ft2m(0, 0),
    ft2m(0.3, 0),
    ft2m(-0.5, 4),
    ft2m(0.2, 0),
]
radii = [
    ft2m(0.1, 0),
    ft2m(1.5, 0),
    ft2m(1.7, 0),
    ft2m(2.7, 0),
    ft2m(2.3, 0),
    ft2m(1, 4),
    ft2m(0.7, 0),
]
shapes = [2, 3, 7, 7, 7, 5, 3]

fuse_sections = [
    FuselageSection(radius = radii[i], center = [xc[i], 0, zc[i]], shape = shapes[i])
    for i in eachindex(xc)
]
fuselage = Fuselage(name = "Main Body", sections = fuse_sections)
AeroGeometry.translate!(fuselage, [ft2m(-5, 0), 0, ft2m(-3, 0)])
#Combine into an airplane
# Create initial Airplane object
airplane = Airplane(
    name = "Cessna 152",
    wings = [wing, horizontal_stabilizer, vertical_stabilizer],
    fuselages = [fuselage],
)

#fig, ax, v = viz(airplane)



plt.plot(airplane, isometric = false)

deflections = Dict(" Aileron " => 20.0, "Elevator" => 20.0, "Rudder" => 20)

# the deflections to the airplane
#deflect_control_surface!(airplane, deflections)

#plt.plot(airplane,isometric=false)
