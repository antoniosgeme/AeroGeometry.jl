# Quickstart

This guide shows you how to build a complete airplane geometry using AeroGeometry.jl. We'll recreate a Cessna 152 by defining its wings, stabilizers, and fuselage.

## Installation

```julia
using Pkg
Pkg.add("AeroGeometry")
```

## Basic Airfoil Operations

First, let's see how to work with airfoils:

```julia
using AeroGeometry
import Plots as plt
plt.plotly()

# Load airfoils from the built-in database
airfoil1 = Airfoil("e221")
airfoil2 = Airfoil("naca6409")

# Blend two airfoils
airfoil3 = blend_airfoils(airfoil1, airfoil2)

# Plot them
plt.plot(airfoil1)
plt.plot!(airfoil2)
plt.plot!(airfoil3)

# Plot the local camber
xc = 0:0.01:1
plt.plot!(xc, local_camber(airfoil3, x_over_c=xc), lw=3)
```

## Building an Airplane

Now let's build a complete airplane. We'll define a helper function to convert feet and inches to meters:

```julia
ft(feet, inches) = 0.3048 * feet + 0.0254 * inches
```

### Define the Main Wing

The main wing is defined by creating wing sections at different spanwise locations:

```julia
wing_sections = [
    WingSection(
        airfoil=Airfoil("naca2412"),
        le_loc=[0, 0, 0],  # Leading edge location [x, y, z]
        chord=ft(5, 4),     # Chord length
    ),
    WingSection(
        airfoil=Airfoil("naca2412"),
        le_loc=[0, ft(7, 0), ft(7, 0) * sind(1)],  # Add dihedral
        chord=ft(5, 4),
    ),
    WingSection(
        airfoil=Airfoil("naca0012"),
        le_loc=[ft(4, 3/4) - ft(3, 8 + 1/2), ft(33, 4)/2, ft(33, 4)/2 * sind(1)],
        chord=ft(3, 8 + 1/2),
        twist=0  # Twist angle in degrees
    )
]

wing = Wing(name="Main Wing", sections=wing_sections, symmetric=true)
```

### Define the Horizontal Stabilizer

```julia
hs_sections = [
    WingSection(
        airfoil=Airfoil("naca0012"),
        le_loc=[0, 0, 0],
        chord=ft(3, 8),
        twist=-2  # Negative twist for tail down force
    ),
    WingSection(
        airfoil=Airfoil("naca0012"),
        le_loc=[ft(1, 0), ft(10, 0) / 2, 0],
        chord=ft(2, 4 + 3 / 8),
        twist=-2
    )
]

horizontal_stabilizer = Wing(
    name="Horizontal Stabilizer", 
    sections=hs_sections, 
    symmetric=true
)

# Move it to the correct position on the aircraft
translate!(horizontal_stabilizer, [4.0648, 0, -0.6096])
```

### Define the Vertical Stabilizer

```julia
vs_sections = [
    WingSection(
        airfoil=Airfoil("naca0012"),
        le_loc=[ft(-5, 0), 0, 0],
        chord=ft(8, 8),
        twist=0
    ),
    WingSection(
        airfoil=Airfoil("naca0012"),
        le_loc=[ft(0, 0), 0, ft(1, 0)],
        chord=ft(3, 8),
        twist=0
    ),
    WingSection(
        airfoil=Airfoil("naca0012"),
        le_loc=[ft(0, 8), 0, ft(5, 0)],
        chord=ft(2, 8), 
        twist=0
    )
]

vertical_stabilizer = Wing(name="Vertical Stabilizer", sections=vs_sections)
translate!(vertical_stabilizer, [ft(16, 11) - ft(3, 8), 0, ft(-2, 0)])
```

### Define the Fuselage

The fuselage is created by defining cross-sections at different longitudinal positions:

```julia
# Define fuselage cross-section positions and dimensions
xc = [0, 0, ft(3, 0), ft(5, 0), ft(10, 4), ft(12, 4), ft(21, 11)]
zc = [ft(-1, 0), ft(-1, 0), ft(-0.85, 0), ft(0, 0), ft(0.3, 0), ft(-0.5, 4), ft(0.2, 0)]
radii = [ft(0.1, 0), ft(1.5, 0), ft(1.7, 0), ft(2.7, 0), ft(2.3, 0), ft(1, 4), ft(0.7, 0)]
shapes = [2, 3, 7, 7, 7, 5, 3]  # Shape parameter for cross-section

fuse_sections = [
    FuselageXSec(radius=radii[i], xyz_c=[xc[i], 0, zc[i]], shape=shapes[i]) 
    for i in eachindex(xc)
]

fuselage = Fuselage(name="Main Body", sections=fuse_sections)
translate!(fuselage, [ft(-5, 0), 0, ft(-3, 0)])
```

### Assemble the Complete Airplane

Finally, combine all components into a single airplane:

```julia
airplane = Airplane(
    name="Cessna 152",
    wings=[wing, horizontal_stabilizer, vertical_stabilizer],
    fuselages=[fuselage]
)

# Visualize it
plt.plot(airplane)
```

![Cessna 152](assets/Cessna152.png)

## Next Steps

Now that you have a complete airplane geometry, you can:
- Export it for use in analysis tools
- Modify dimensions and airfoil shapes
- Add control surfaces
- Create mesh representations

Check out the [Reference](@ref reference) documentation for details on all available functions and types.
