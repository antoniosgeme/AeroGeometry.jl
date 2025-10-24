# AeroGeometry.jl

Documentation for AeroGeometry.jl - A Julia package for aerodynamic geometry modeling.

## Features

- Airfoil geometry and operations
- Wing and fuselage modeling
- Airplane assembly
- Visualization with Plots.jl and Makie.jl

## Installation

```julia
using Pkg
Pkg.add("AeroGeometry")
```

## Quick Start

```julia
using AeroGeometry

# Load an airfoil from the database
airfoil = Airfoil("naca0012")

# Create a wing
wing = Wing([
    WingSection(airfoil, [0.0, 0.0, 0.0], 1.0, 0.0)
])
```
