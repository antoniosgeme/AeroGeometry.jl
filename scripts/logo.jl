using AeroGeometry
using Plots
plotly()
using Colors

JULIA_PURPLE = RGB(0.584, 0.345, 0.698)
JULIA_RED = RGB(0.796, 0.235, 0.200)
JULIA_GREEN = RGB(0.220, 0.596, 0.149)
JULIA_BLUE = RGB(0.251, 0.388, 0.847)
cessna = cessna152()

# Manual plotting with full color control
p = plot(legend=false, ticks=false, border=:none, 
         xlabel="", ylabel="", zlabel="",
         camera=(45, 30));

# Plot Main Wing (wing 1) in PURPLE
wing1 = deepcopy(cessna.wings[1])
x1, y1, z1 = AeroGeometry.coordinates(wing1)
# Plot wing sections
for i in axes(x1, 2)
    plot!(p, x1[:, i], y1[:, i], z1[:, i], color=JULIA_BLUE, lw=3);
    if wing1.symmetric && i > 1
        plot!(p, x1[:, i], y1[:, 1] .- y1[:, i], z1[:, i], color=JULIA_BLUE, lw=3);
    end
end
# Plot surface
plot!(p, x1, y1, z1, seriestype=:surface, color=JULIA_PURPLE, alpha=1);
if wing1.symmetric
    plot!(p, x1, y1[:, 1] .- y1, z1, seriestype=:surface, color=JULIA_PURPLE, alpha=1);
end

# Plot Horizontal Stabilizer (wing 2) in GREEN
wing2 = deepcopy(cessna.wings[2])
x2, y2, z2 = AeroGeometry.coordinates(wing2)
# Plot wing sections
for i in axes(x2, 2)
    plot!(p, x2[:, i], y2[:, i], z2[:, i], color=JULIA_BLUE, lw=3);
    if wing2.symmetric && i > 1
        plot!(p, x2[:, i], y2[:, 1] .- y2[:, i], z2[:, i], color=JULIA_BLUE, lw=3);
    end
end
# Plot surface
plot!(p, x2, y2, z2, seriestype=:surface, color=JULIA_GREEN, alpha=1);
if wing2.symmetric
    plot!(p, x2, y2[:, 1] .- y2, z2, seriestype=:surface, color=JULIA_GREEN, alpha=1)
end

# Plot Vertical Stabilizer (wing 3) in BLUE
wing3 = deepcopy(cessna.wings[3])
x3, y3, z3 = AeroGeometry.coordinates(wing3)
# Plot wing sections
for i in axes(x3, 2)
    plot!(p, x3[:, i], y3[:, i], z3[:, i], color=JULIA_BLUE, lw=3)
    if wing3.symmetric && i > 1
        plot!(p, x3[:, i], y3[:, 1] .- y3[:, i], z3[:, i], color=JULIA_BLUE, lw=3)
    end
end
# Plot surface
plot!(p, x3, y3, z3, seriestype=:surface, color=JULIA_GREEN, alpha=1)
if wing3.symmetric
    plot!(p, x3, y3[:, 1] .- y3, z3, seriestype=:surface, color=JULIA_GREEN, alpha=1)
end

# Plot Fuselage in RED
fuselage = cessna.fuselages[1]
# Plot fuselage cross-sections
for xsec in fuselage.sections
    x, y, z = AeroGeometry.coordinates(xsec)
    push!(x, first(x))  # Close the loop
    push!(y, first(y))
    push!(z, first(z))
    plot!(p, x, y, z, color=JULIA_BLUE, lw=3)
end
# Plot fuselage surface
if length(fuselage.sections) > 1
    num_points = length(AeroGeometry.coordinates(fuselage.sections[1])[1])
    num_sections = length(fuselage.sections)
    
    x_fuse = zeros(num_points+1, num_sections)
    y_fuse = zeros(num_points+1, num_sections)
    z_fuse = zeros(num_points+1, num_sections)
    
    for i = 1:num_sections
        x_fuse[1:(end-1), i], y_fuse[1:(end-1), i], z_fuse[1:(end-1), i] =
            AeroGeometry.coordinates(fuselage.sections[i])
        x_fuse[end, i] = x_fuse[1, i]  # Close the loop
        y_fuse[end, i] = y_fuse[1, i]
        z_fuse[end, i] = z_fuse[1, i]
    end
    
    plot!(p, x_fuse, y_fuse, z_fuse, seriestype=:surface, color=JULIA_RED, alpha=1)
end

display(p)