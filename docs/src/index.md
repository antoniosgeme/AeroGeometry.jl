```@meta
CurrentModule = AeroGeometry
```

# AeroGeometry

Documentation for [AeroGeometry](https://github.com/antoniosgeme/AeroGeometry.jl).

A helper library to import and manipulate aircraft geometry.

## Installation

```julia
using Pkg
Pkg.add("https://github.com/antoniosgeme/AeroGeometry.jl.git")
```

## Quick Start

### Working with Airfoils

Load and manipulate airfoils from the built-in database:

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

camb = camber(airfoil3)
thick = thickness(airfoil3)
```

### Building Aircraft Geometry

Create wings by defining sections at different spanwise locations:

```julia
using AeroGeometry: ft2m

wing_sections = [
    WingSection(
        airfoil=Airfoil("naca2412"),
        le_loc=[0, 0, 0],  # Leading edge location [x, y, z]
        chord=ft2m(5, 4),  # Chord length in meters
    ),
    WingSection(
        airfoil=Airfoil("naca0012"),
        le_loc=[ft2m(4, 3/4) - ft2m(3, 8 + 1/2), ft2m(33, 4)/2, ft2m(33, 4)/2 * sind(1)],
        chord=ft2m(3, 8 + 1/2),
        twist=0  # Twist angle in degrees
    )
]

wing = Wing(name="Main Wing", sections=wing_sections, symmetric=true)
plt.plot(wing)
```

### Complete Airplane

Combine multiple components into a complete airplane:

```julia
airplane = Airplane(
    name="Cessna 152",
    wings=[wing, horizontal_stabilizer, vertical_stabilizer],
    fuselages=[fuselage]
)

plt.plot(airplane)
```

## Features

- Import airfoils from UIUC database or custom files
- Generate NACA 4-digit airfoils
- Blend and manipulate airfoil geometries
- Define wings with variable sweep, dihedral, twist, and taper
- Create fuselages with custom cross-sections
- Assemble complete aircraft configurations
- Visualization with Plots.jl
- Export geometries for analysis tools

