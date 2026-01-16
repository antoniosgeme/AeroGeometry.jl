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

function boeing737()
    # Dimensions based on Boeing 737-800 specs (length 39.47 m, span 35.79 m,
    # fuselage width 3.76 m, wing area 124.6 m^2). Tail sizes are approximate.
    length = 39.47
    wing_span = 35.79
    wing_area = 124.6
    fuselage_diameter = 3.76

    # Main wing
    wing_root_chord = 5.5
    wing_tip_chord = 2 * wing_area / wing_span - wing_root_chord
    wing_sweep = 25.0
    wing_dihedral = 6.0
    wing_le_x = 12.5
    wing_z = -0.6

    wing_sections = [
        WingSection(
            airfoil = Airfoil("naca2412"),
            position = [y * tand(wing_sweep), y, y * tand(wing_dihedral)],
            chord = wing_root_chord -
                (wing_root_chord - wing_tip_chord) * (y / (wing_span / 2)),
        )
        for y in range(0, stop = wing_span / 2, length = 6)
    ]
    wing = Wing(name = "Main Wing", sections = wing_sections, symmetric = true)
    translate!(wing, [wing_le_x, 0, wing_z])

    # Horizontal stabilizer (approximate)
    hs_span = 12.6
    hs_area = 32.7
    hs_root_chord = 3.2
    hs_tip_chord = 2 * hs_area / hs_span - hs_root_chord
    hs_sweep = 30.0
    hs_dihedral = 5.0
    hs_le_x = 31.5
    hs_z = 0.8

    hs_sections = [
        WingSection(
            airfoil = Airfoil("naca0012"),
            position = [y * tand(hs_sweep), y, y * tand(hs_dihedral)],
            chord = hs_root_chord -
                (hs_root_chord - hs_tip_chord) * (y / (hs_span / 2)),
            twist = 0,
        )
        for y in range(0, stop = hs_span / 2, length = 4)
    ]
    horizontal_stabilizer =
        Wing(name = "Horizontal Stabilizer", sections = hs_sections, symmetric = true)
    translate!(horizontal_stabilizer, [hs_le_x, 0, hs_z])

    # Vertical stabilizer (approximate)
    vs_span = 6.5
    vs_area = 26.0
    vs_root_chord = 5.5
    vs_tip_chord = 2 * vs_area / vs_span - vs_root_chord
    vs_sweep = 35.0
    vs_le_x = 30.5
    vs_z = fuselage_diameter / 2

    vs_sections = [
        WingSection(
            airfoil = Airfoil("naca0012"),
            position = [z * tand(vs_sweep), 0, z],
            chord = vs_root_chord -
                (vs_root_chord - vs_tip_chord) * (z / vs_span),
            twist = 0,
        )
        for z in range(0, stop = vs_span, length = 4)
    ]
    vertical_stabilizer =
        Wing(name = "Vertical Stabilizer", sections = vs_sections, symmetric = false)
    translate!(vertical_stabilizer, [vs_le_x, 0, vs_z])

    # Fuselage (approximate)
    fuse_radius = fuselage_diameter / 2
    xc = [0.0, 1.5, 4.0, 8.0, 20.0, 30.0, 36.0, length]
    zc = [0.0 for _ in xc]
    radii = [0.3, 1.2, 1.7, fuse_radius, fuse_radius, 1.5, 0.7, 0.2]
    exponents = [2, 3, 4, 4, 4, 4, 3, 2]

    fuse_sections = [
        FuselageSection(
            center = [xc[i], 0, zc[i]],
            shape = Superellipse(2 * radii[i], 2 * radii[i], exponents[i]),
        )
        for i in eachindex(xc)
    ]
    fuselage = Fuselage(name = "Main Body", sections = fuse_sections)

    airplane = Airplane(
        name = "Boeing 737-800",
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

