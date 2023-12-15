# AeroGeometry.jl
A helper library to import and manipulate airfoil and wing geometry.

```julia
using AeroGeometry
using Plots
airfoil1 = Airfoil("naca0012")
airfoil2 = Airfoil("naca6409")
airfoil3 = blend_airfoils(airfoil1,airfoil2)
scatter(airfoil)
scatter!(airfoil2)
scatter!(airfoil3)
