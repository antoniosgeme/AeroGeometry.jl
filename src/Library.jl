function cessna152()

    # Define the wings
    wing_sections = [
        WingSection(airfoil = Airfoil("naca2412"), position = [0, 0, 0], chord = ft2m(5, 4)),
        WingSection(
            airfoil = Airfoil("naca2412"),
            position = [0, ft2m(7), ft2m(7) * sind(1)],
            chord = ft2m(5, 4),
        ),
        WingSection(
            airfoil = Airfoil("naca0012"),
            position = [
                ft2m(4, 3/4) - ft2m(3, 8 + 1/2),
                ft2m(33, 4)/2,
                ft2m(33, 4)/2 * sind(1),
            ],
            chord = ft2m(3, 8 + 1/2),
            twist = 0,
        ),
    ]
    wing = Wing(name = "Main Wing", sections = wing_sections, symmetric = true)

    hs_sections = [
        WingSection(
            airfoil = Airfoil("naca0012"),
            position = [0, 0, 0],
            chord = ft2m(3, 8),
            twist = -2,
        ),
        WingSection(
            airfoil = Airfoil("naca0012"),
            position = [ft2m(1), ft2m(10) / 2, 0],
            chord = ft2m(2, 4 + 3 / 8),
            twist = -2,
        ),
    ]
    horizontal_stabilizer =
        Wing(name = "Horizontal Stabilizer", sections = hs_sections, symmetric = true)
    translate!(horizontal_stabilizer, [4.0648, 0, -0.6096])

    vs_sections = [
        WingSection(
            airfoil = Airfoil("naca0012"),
            position = [ft2m(-5), 0, 0],
            chord = ft2m(8, 8),
            twist = 0,
        ),
        WingSection(
            airfoil = Airfoil("naca0012"),
            position = [0, 0, ft2m(1)],
            chord = ft2m(3, 8),
            twist = 0,
        ),
        WingSection(
            airfoil = Airfoil("naca0012"),
            position = [ft2m(0, 8), 0, ft2m(5)],
            chord = ft2m(2, 8),
            twist = 0,
        ),
    ]
    vertical_stabilizer =
        Wing(name = "Vertical Stabilizer", sections = vs_sections, symmetric = false)
    translate!(vertical_stabilizer, [ft2m(16, 11) - ft2m(3, 8), 0, ft2m(-2)])

    # Define the fuselage
    xc = [0, 0, ft2m(3), ft2m(5), ft2m(10, 4), ft2m(12, 4), ft2m(21, 11)]
    zc = [ft2m(-1), ft2m(-1), ft2m(-0.85), ft2m(0), ft2m(0.3), ft2m(-0.5, 4), ft2m(0.2)]
    radii = [ft2m(0.1), ft2m(1.5), ft2m(1.7), ft2m(2.7), ft2m(2.3), ft2m(1, 4), ft2m(0.7)]
    exponents = [2, 3, 7, 7, 7, 5, 3]

    fuse_sections = [
        FuselageSection(
            center = [xc[i], 0, zc[i]], 
            shape = Superellipse(2 * radii[i], 2 * radii[i], exponents[i])
        )
        for i in eachindex(xc)
    ]
    fuselage = Fuselage(name = "Main Body", sections = fuse_sections)
    translate!(fuselage, [ft2m(-5), 0, ft2m(-3)])
    # Combine into an airplane
    airplane = Airplane(
        name = "Cessna 152",
        wings = [wing, horizontal_stabilizer, vertical_stabilizer],
        fuselages = [fuselage],
    )

    return airplane
end

function boeing777()
    # Boeing 777 geometry based on AeroFuse tutorials-stability.md
    fuselage_radius = 3.04
    fuselage_length = 63.7
    x_a = 0.15
    x_b = 0.7
    c_nose = 2.0
    c_rear = 1.2
    d_nose = -0.5
    d_rear = 1.0

    # Main wing (two-span segments, three chord stations)
    wing_chords = [14.0, 9.73, 1.43561]
    wing_spans = [14.0, 46.9] ./ 2
    wing_dihedral = 6.0
    wing_sweep = 35.6
    wing_position = [19.51, 0.0, -2.5]
    wing_airfoils = [Airfoil("rae2822"), Airfoil("rae2822"), Airfoil("naca0012")]

    wing_y = cumsum(vcat(0.0, wing_spans))
    wing_sections = [
        WingSection(
            airfoil = wing_airfoils[i],
            position = [wing_y[i] * tand(wing_sweep), wing_y[i], wing_y[i] * tand(wing_dihedral)],
            chord = wing_chords[i],
        )
        for i in eachindex(wing_y)
    ]
    wing = Wing(name = "Main Wing", sections = wing_sections, symmetric = true)
    translate!(wing, wing_position)

    # Tail sizing helper
    function trapezoid_planform(area::Real, aspect::Real, taper::Real)
        span = sqrt(area * aspect)
        root_chord = 2 * area / (span * (1 + taper))
        tip_chord = root_chord * taper
        return span, root_chord, tip_chord
    end

    # Horizontal stabilizer
    hs_area = 101.0
    hs_aspect = 4.2
    hs_taper = 0.4
    hs_dihedral = 7.0
    hs_sweep = 35.0
    hs_span, hs_root_chord, hs_tip_chord = trapezoid_planform(hs_area, hs_aspect, hs_taper)
    hs_y = [0.0, hs_span / 2]
    hs_sections = [
        WingSection(
            airfoil = Airfoil("naca0012"),
            position = [hs_y[i] * tand(hs_sweep), hs_y[i], hs_y[i] * tand(hs_dihedral)],
            chord = i == 1 ? hs_root_chord : hs_tip_chord,
        )
        for i in eachindex(hs_y)
    ]
    horizontal_stabilizer =
        Wing(name = "Horizontal Stabilizer", sections = hs_sections, symmetric = true)
    htail_position = [fuselage_length - 8.0, 0.0, 0.0]
    translate!(horizontal_stabilizer, htail_position)
    rotate!(horizontal_stabilizer, axis = [0, 1, 0], angle = -2.0)

    # Vertical stabilizer
    vs_area = 56.1
    vs_aspect = 1.5
    vs_taper = 0.4
    vs_sweep = 44.4
    vs_span, vs_root_chord, vs_tip_chord = trapezoid_planform(vs_area, vs_aspect, vs_taper)
    vs_z = [0.0, vs_span]
    vs_sections = [
        WingSection(
            airfoil = Airfoil("naca0009"),
            position = [vs_z[i] * tand(vs_sweep), 0.0, vs_z[i]],
            chord = i == 1 ? vs_root_chord : vs_tip_chord,
        )
        for i in eachindex(vs_z)
    ]
    vertical_stabilizer =
        Wing(name = "Vertical Stabilizer", sections = vs_sections, symmetric = false)
    vtail_position = htail_position .- [2.0, 0.0, 0.0]
    translate!(vertical_stabilizer, vtail_position)

    # Fuselage (hyperellipse-based, sampled into sections)
    hyperellipse(ξ, a) = (1 - ξ^a)^(1 / a)
    n_nose = 6
    n_cabin = 4
    n_rear = 6
    ts_nose = range(0, stop = 1, length = n_nose)
    ts_cabin = range(0, stop = 1, length = n_cabin)
    ts_rear = range(0, stop = 1, length = n_rear)

    L_nose = x_a * fuselage_length
    L_cabin = (x_b - x_a) * fuselage_length
    L_rear = (1 - x_b) * fuselage_length

    x_nose = ts_nose .* L_nose
    x_cabin = ts_cabin .* L_cabin .+ x_nose[end]
    x_rear = ts_rear .* L_rear .+ x_cabin[end]

    r_nose = fuselage_radius .* reverse(hyperellipse.(ts_nose, c_nose))
    r_cabin = fill(fuselage_radius, length(ts_cabin))
    r_rear = fuselage_radius .* hyperellipse.(ts_rear, c_rear)

    z_nose = collect(range(d_nose, stop = 0.0, length = length(x_nose)))
    z_cabin = fill(0.0, length(x_cabin))
    z_rear = collect(range(0.0, stop = d_rear, length = length(x_rear)))

    xc = vcat(x_nose, x_cabin[2:end], x_rear[2:end])
    radii = vcat(r_nose, r_cabin[2:end], r_rear[2:end])
    zc = vcat(z_nose, z_cabin[2:end], z_rear[2:end])

    fuse_sections = [
        FuselageSection(
            center = [xc[i], 0.0, zc[i]],
            shape = Superellipse(2 * radii[i]+0.01, 2 * radii[i]+0.01, 2.0),
        )
        for i in eachindex(xc)
    ]
    fuselage = Fuselage(name = "Main Body", sections = fuse_sections)

    airplane = Airplane(
        name = "Boeing 777",
        wings = [wing, horizontal_stabilizer, vertical_stabilizer],
        fuselages = [fuselage],
    )

    return airplane
end



function rectangular_wing(
    span::Real,
    chord::Real;
    airfoil=Airfoil("naca2412"),
    name::String = "Rectangular Wing",
    n_sections::Int = 5,
    symmetric::Bool = true,
)
    section_positions = range(0, stop = span / 2, length = n_sections)
    sections = [
        WingSection(airfoil = airfoil, position = [0, y, 0], chord = chord)
        for y in section_positions
    ]
    wing = Wing(name = name, sections = sections, symmetric = symmetric)
    return wing
end

function tapered_wing(
    span::Real,
    root_chord::Real,
    tip_chord::Real;
    airfoil=Airfoil("naca2412"),
    name::String = "Tapered Wing",
    n_sections::Int = 5,
    symmetric::Bool = true,
    sweep::Real = 0.0,  # Sweep angle in degrees
)
    section_positions = range(0, stop = span / 2, length = n_sections)
    
    # Calculate sweep offset (assumes sweep is given as distance at tip)
    # If you want to use sweep angle instead, use: x_offset = y * tand(sweep)
    sections = [
        WingSection(
            airfoil = airfoil,
            position = [y * tand(sweep), y, 0],  # sweep in degrees,  # Linear sweep
            chord = root_chord - (root_chord - tip_chord) * (y / (span / 2)),
        )
        for y in section_positions
    ]
    wing = Wing(name = name, sections = sections, symmetric = symmetric)
    return wing
end

