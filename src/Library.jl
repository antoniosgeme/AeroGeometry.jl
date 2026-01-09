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
    shapes = [2, 3, 7, 7, 7, 5, 3]

    fuse_sections = [
        FuselageSection(radius = radii[i], center = [xc[i], 0, zc[i]], shape = shapes[i])
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