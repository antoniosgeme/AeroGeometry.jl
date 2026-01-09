rectangular_wing() = Wing(
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
    wing = rectangular_wing()

    @test area(wing) ≈ 4
    @test area(wing,centerline=true) ≈ 10

end