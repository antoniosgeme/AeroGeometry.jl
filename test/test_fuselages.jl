@testset "Fuselage shapes" begin
    shape1 = Superellipse(2.0, 2.0, 2.0)
    shape2 = Superellipse(4.0, 4.0, 2.0)
    @test area(shape2) ≈ 4 * area(shape1) atol=1e-6

    square = [
        -1.0 -1.0;
         1.0 -1.0;
         1.0  1.0;
        -1.0  1.0;
    ]
    section = FuselageSection(shape = square)
    @test section.shape isa ArbitraryShape
    @test area(section.shape) ≈ 4.0 atol=1e-6
end

@testset "Fuselage geometry" begin
    section1 = FuselageSection(center = [0.0, 0.0, 0.0], shape = Superellipse(2.0, 2.0, 2.0))
    section2 = FuselageSection(center = [10.0, 0.0, 0.0], shape = Superellipse(2.0, 2.0, 2.0))
    fuselage = Fuselage(name = "Test Fuselage", sections = [section1, section2])

    @test length(fuselage) ≈ 10.0 atol=1e-6
    @test volume(fuselage) ≈ 10.0 * area(section1) atol=1e-6

    θ = 0:0.1:2π
    X, Y, Z = coordinates(fuselage; θ = θ)
    @test size(X) == (length(θ), 2)
    @test size(Y) == (length(θ), 2)
    @test size(Z) == (length(θ), 2)

    @test area(fuselage, closed = true) > area(fuselage)
end

@testset "Fuselage transforms" begin
    fuselage = Fuselage(sections = [FuselageSection(), FuselageSection(center = [2.0, 0.0, 0.0])])
    translate!(fuselage, [1.0, 0.0, 0.0])
    @test fuselage.position == [1.0, 0.0, 0.0]

    rotate!(fuselage, axis = [0.0, 0.0, 1.0], angle = 90.0, pivot = [-1.0, 0.0, 0.0])
    @test isapprox(fuselage.position[1], 0.0; atol = 1e-6)
    @test isapprox(fuselage.position[2], 1.0; atol = 1e-6)
end
