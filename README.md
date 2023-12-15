# AeroGeometry.jl
A helper library to import and manipulate airfoil and wing geometry.

```julia
using AeroGeometry
using Plots

airfoil1 = Airfoil("e221")
airfoil2 = Airfoil("naca6409")
airfoil3 = blend_airfoils(airfoil1,airfoil2)

plot(airfoil1)
plot!(airfoil2)
plot!(airfoil3)

xc = 0:0.01:1
plot!(xc,get_local_camber(airfoil3,x_over_c=xc),lw=3)


