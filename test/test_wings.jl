test_wing() = Wing(
    name = "Rectangular Wing",
    sections = [
        WingSection(
            position = [0.0, 3, 0.0],
            chord = 1.0,
            airfoil = Airfoil("NACA0012"),
            twist = 0.0,
        ),
        WingSection(
            position = [0.0, 5.0, 0.0],
            chord= 1.0,
            airfoil = Airfoil("NACA0012"),
            twist = 0.0,
        )
    ],
    symmetric = true
)

@testset "Area function tests" begin
    wing = test_wing()

    @test area(wing) ≈ 4
    @test area(wing,centerline=true) ≈ 10


    wing = Wing()
    @test area(wing) ≈ 2
    @test area(wing,centerline=true) ≈ 2
end

function tapered_test_wing(; root_y = 0.0)
    return Wing(
        name = "Tapered Wing",
        sections = [
            WingSection(
                position = [0.0, root_y, 0.0],
                chord = 2.0,
                airfoil = Airfoil("NACA0012"),
                twist = 0.0,
            ),
            WingSection(
                position = [0.0, root_y + 5.0, 0.0],
                chord = 1.0,
                airfoil = Airfoil("NACA0012"),
                twist = 0.0,
            ),
        ],
        symmetric = true,
    )
end

@testset "Wing geometry helpers" begin
    wing = tapered_test_wing()

    expected_span = 2 * norm(quarter_chord(wing.sections[2]) - quarter_chord(wing.sections[1]))
    @test span(wing) ≈ expected_span atol=1e-6
    @test area(wing) ≈ 15.0 atol=1e-1
    @test aspect_ratio(wing) ≈ (expected_span^2 / 15.0) atol=1e-2
    @test taper_ratio(wing) ≈ 0.5 atol=1e-6
    @test mean_geometric_chord(wing) ≈ 1.5 atol=1e-2
    @test mean_aerodynamic_chord(wing) ≈ 1.5556 atol=1e-3

    ac = aerodynamic_center(wing)
    @test abs(ac[2]) < 1e-8
    @test ac[1] > 0

    X, Y, Z = coordinates(wing, camberline = true)
    @test size(X, 2) == 3
    @test size(Y, 2) == 3
    @test size(Z, 2) == 3
end

@testset "Wing span with centerline" begin
    wing = tapered_test_wing(root_y = 1.0)
    expected_span = 2 * norm(quarter_chord(wing.sections[2]) - quarter_chord(wing.sections[1]))
    centerline_span = expected_span + 2 * abs(quarter_chord(wing.sections[1])[2])
    @test span(wing) ≈ expected_span atol=1e-6
    @test span(wing, centerline = true) ≈ centerline_span atol=1e-6
end

@testset "Wing transforms and validation" begin
    wing = tapered_test_wing()
    translate!(wing, [1.0, 2.0, 3.0])
    @test wing.position == [1.0, 2.0, 3.0]

    rotate!(
        wing,
        axis = [0.0, 0.0, 1.0],
        angle = 90.0,
        pivot = [-1.0, -2.0, -3.0],
    )
    @test isapprox(wing.orientation[1, 1], 0.0; atol = 1e-6)

    bad_wing = Wing(
        sections = [
            WingSection(position = [0.0, -1.0, 0.0], chord = 1.0, airfoil = Airfoil("NACA0012")),
            WingSection(position = [0.0, 2.0, 0.0], chord = 1.0, airfoil = Airfoil("NACA0012")),
        ],
        symmetric = false,
    )
    bad_wing.symmetric = true
    @test_throws ErrorException AeroGeometry.validate!(bad_wing)
end
