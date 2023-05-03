using AeroFlux
using LinearAlgebra
using Plots
# Options
aoa = 1 # Angle of attack in radians
U_inf = 1 # Total freestream velocity
N = 100 #Number of points on airfoil surfaces (2N-1 total)
airfoil = Airfoil("naca6412")

x_vort = airfoil.coordinates[:,1]
y_vort = airfoil.coordinates[:,2]
x_col = ( x_vort[2:end] + x_vort[1:end-1] ) /2
y_col = ( y_vort[2:end] + y_vort[1:end-1] ) /2
panel_lengths = hypot.(diff(x_vort),diff(y_vort))
panel_angles = atan.(diff(y_vort),diff(x_vort))
panel_normals = [-sin.(panel_angles) cos.(panel_angles)]

TE_length = hypot(x_vort[1] - x_vort[end],y_vort[1]-y_vort[end])
TE_angle = atan(y_vort[1]-y_vort[end], x_vort[1] - x_vort[end])
upper_TE_vec = airfoil.coordinates[1,:]
lower_TE_vec = airfoil.coordinates[end,:]
tangent_TE_vec = upper_TE_vec - lower_TE_vec
tangent_TE_vec = tangent_TE_vec / hypot(tangent_TE_vec[1],tangent_TE_vec[2])
end_panel_tangent1 = [ -cos( panel_angles[1] ), -sin( panel_angles[1] ) ]
end_panel_tangent2 = [ cos( panel_angles[end] ), sin( panel_angles[end] ) ]
tangent_angle1 = atan(end_panel_tangent1[2],end_panel_tangent1[1])
tangent_angle2 = atan(end_panel_tangent2[2],end_panel_tangent2[1])
bisect_TE_vec = [ real(exp(1im*(tangent_angle1+tangent_angle2)/2)),imag( exp(1im*(tangent_angle1+tangent_angle2)/2))]
cross_prod = abs(bisect_TE_vec[2]*tangent_TE_vec[1] - bisect_TE_vec[1]*tangent_TE_vec[2])
dot_prod = abs(dot(bisect_TE_vec,tangent_TE_vec))

fig = plotme(airfoil)
quiver!(fig,x_vort,y_vort,quiver=(panel_normals[:,1],panel_normals[:,2]))


