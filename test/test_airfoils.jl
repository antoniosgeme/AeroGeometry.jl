using Test
include("../src/Tools.jl")  # Adjust the path if needed
include("../src/Airfoils.jl")  # Adjust the path if needed

function create_test_airfoil()
    coords = [
        1.0000  0.0013;
        0.9750  0.0047;
        0.9026  0.0141;
        0.7900  0.0273;
        0.6485  0.0415;
        0.4921  0.0534;
        0.3365  0.0597;
        0.1972  0.0572;
        0.0882  0.0447;
        0.0203  0.0237;
        0.0003 -0.0028;
        0.0302 -0.0285;
        0.1070 -0.0480;
        0.2230 -0.0585;
        0.3668 -0.0591;
        0.5238 -0.0514;
        0.6784 -0.0387;
        0.8153 -0.0245;
        0.9206 -0.0119;
        0.9840 -0.0035;
        1.0000  -0.0013
    ]
    return Airfoil("test_airfoil", coords[:,1],coords[:,2])
end

@testset "Airfoil Initialization" begin
    airfoil = create_test_airfoil()

    @test airfoil.name == "test_airfoil"
    @test size(coordinates(airfoil)) == (21,2)

end

@testset "Coordinate Functions" begin
    airfoil = create_test_airfoil()

    @test leading_edge_index(airfoil) == 11  # Ensure it finds LE at x=0
    @test size(coordinates(airfoil, :all)) == (21, 2)
    @test size(coordinates(airfoil, :upper)) == (11, 2)
    @test size(coordinates(airfoil, :lower)) == (11, 2)

    @test trailing_edge_thickness(airfoil) ≈ 0.0026  # Test simple case
end

@testset "Airfoil Geometry" begin
    airfoil = create_test_airfoil()

    @test max_camber(airfoil) > 0
    @test max_thickness(airfoil) > 0
    @test area(airfoil) > 0

    centroid_x, centroid_y = centroid(airfoil)
    @test centroid_x > 0
end



@testset "Repaneling" begin
    airfoil = create_test_airfoil()

    repaneled = repanel(airfoil, 10)
    @test length(repaneled.x) == 19
end

@testset "Airfoil Blending" begin
    airfoil1 = create_test_airfoil()
    airfoil2 = Airfoil("NACA6409")

    blended = blend_airfoils(airfoil1, airfoil2, fraction=0.5, points_per_side=10)
    @test blended.name == "test_airfoil+NACA6409"
    @test size(coordinates(blended)) == (19, 2)
end
