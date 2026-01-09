
@testset "README examples" begin
    @testset "Airfoil blending" begin
        airfoil1 = Airfoil("naca2412")
        airfoil2 = Airfoil("e221")
        blended = blend_airfoils(airfoil1, airfoil2, fraction=0.5)
        @test blended isa Airfoil
        camb = camber(blended)
        @test size(camb, 2) == 2  # Should have x,y columns
    end
    
    @testset "Wing creation" begin
        wing = Wing(
            name="Main Wing",
            sections=[
                WingSection(airfoil=Airfoil("naca2412"), position=[0, 0, 0], chord=1.524),
                WingSection(airfoil=Airfoil("naca0012"), position=[1.2, 5.0, 0.15], chord=1.1, twist=-2.0)
            ],
            symmetric=true
        )
        @test wing isa Wing
        @test length(wing.sections) == 2
        @test wing.symmetric == true
    end
end