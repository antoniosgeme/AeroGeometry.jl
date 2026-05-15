

function create_test_airfoil()
    coords = [
        1.0000 0.0013;
        0.9750 0.0047;
        0.9026 0.0141;
        0.7900 0.0273;
        0.6485 0.0415;
        0.4921 0.0534;
        0.3365 0.0597;
        0.1972 0.0572;
        0.0882 0.0447;
        0.0203 0.0237;
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
        1.0000 -0.0013
    ]
    return Airfoil("test_airfoil", coords[:, 1], coords[:, 2])
end

@testset "Airfoil Initialization" begin
    airfoil = create_test_airfoil()

    @test airfoil.name == "test_airfoil"
    @test size(coordinates(airfoil)) == (21, 2)

end

@testset "NACA Constructor Options" begin
    default_te = Airfoil("NACA0012")
    sharp_te = Airfoil("NACA0012"; te_sharp = true)
    lower_resolution = Airfoil("NACA0012"; points_per_side = 40)

    @test trailing_edge_thickness(default_te) > 0
    @test trailing_edge_thickness(sharp_te) ≈ 0 atol = 1e-12
    @test length(lower_resolution.x) == 79
end

@testset "Ignored NACA Options" begin
    @test_logs (:warn, r"te_sharp only applies.*UIUC") Airfoil("e221"; te_sharp = true)

    mktempdir() do dir
        path = joinpath(dir, "test_airfoil.dat")
        open(path, "w") do io
            write(io, "test_airfoil\n")
            write(io, "1.0 0.01\n")
            write(io, "0.0 0.0\n")
            write(io, "1.0 -0.01\n")
        end

        @test_logs (:warn, r"te_sharp only applies.*file-loaded") Airfoil(path; te_sharp = true)
    end
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

    @test max_camber(airfoil) ≈ 1
    @test max_thickness(airfoil) ≈ 1
    @test area(airfoil) ≈ 0.08084719

    centroid_x, centroid_y = centroid(airfoil)
    @test centroid_x ≈ 0.42324898501935143
    @test centroid_y ≈ -1.0646356416345823e-5
end



@testset "Repaneling" begin
    airfoil = create_test_airfoil()

    repaneled = repanel(airfoil, 30)
    @test length(repaneled.x) == 30
end

@testset "Airfoil Blending" begin
    airfoil1 = create_test_airfoil()
    airfoil2 = Airfoil("NACA6409")

    blended = blend_airfoils(airfoil1, airfoil2, fraction = 0.5, points_per_side = 10)
    @test blended.name == "test_airfoil+NACA6409"
end

@testset "Normalization and Utilities" begin
    airfoil = create_test_airfoil()
    airfoil.x .= airfoil.x .* 2 .+ 5
    airfoil.y .= airfoil.y .* 2 .+ 1
    normalize!(airfoil)

    @test minimum(airfoil.x) ≈ 0 atol=1e-5
    @test maximum(airfoil.x) ≈ 1 atol=1e-5
    @test abs((airfoil.y[1] + airfoil.y[end]) / 2) < 1e-6

    tangents_center = tangents(airfoil, centers = true)
    tangents_nodes = tangents(airfoil, centers = false)
    normals_nodes = normals(airfoil, centers = false)

    @test size(tangents_center, 2) == 2
    @test size(tangents_nodes, 2) == 2
    @test size(normals_nodes, 2) == 2
    @test size(tangents_nodes, 1) == length(airfoil.x)
end

@testset "Control Surface Deflection" begin
    airfoil = create_test_airfoil()
    deflected = deflect_control_surface(airfoil, deflection = 10.0, x_hinge = 0.75)

    @test deflected isa Airfoil
    @test sum(abs.(deflected.y .- airfoil.y)) > 0
end

@testset "Airfoil Database Listing" begin
    names = list_airfoil_names("a18")
    @test "a18" in names
end
