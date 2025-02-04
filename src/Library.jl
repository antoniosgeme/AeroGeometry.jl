function Cessna152()

    # Define the wings
    wing_xsecs = [
        WingXSec(
            airfoil=Airfoil("naca2412"),
            le_loc=[0, 0, 0],
            chord=ft2m(5,4),
        ),
        WingXSec(
            airfoil=Airfoil("naca2412"),
            le_loc=[0, ft2m(7), ft2m(7) * sind(1)],
            chord=ft2m(5,4),
        ),
        WingXSec(
            airfoil=Airfoil("naca0012"),
            le_loc=[ft2m(4, 3/4) - ft2m(3, 8 + 1/2), ft2m(33, 4)/2, ft2m(33, 4)/2 * sind(1)],
            chord=ft2m(3, 8 + 1/2),
            twist=0
        )
    ]
    wing = Wing(name="Main Wing", xsecs=wing_xsecs, symmetric=true)

    hs_xsecs = [
        WingXSec(
            airfoil=Airfoil("naca0012"),
            le_loc=[0, 0, 0],
            chord=ft2m(3,8),
            twist=-2
        ),
        WingXSec(
            airfoil=Airfoil("naca0012"),
            le_loc=[ft2m(1), ft2m(10) / 2, 0],
            chord=ft2m(2, 4 + 3 / 8),
            twist=-2
        )
    ]
    horizontal_stabilizer = Wing(name="Horizontal Stabilizer", xsecs=hs_xsecs, symmetric=true)
    translate!(horizontal_stabilizer, [4.0648, 0, -0.6096])

    vs_xsecs = [
        WingXSec(
            airfoil=Airfoil("naca0012"),
            le_loc=[ft2m(-5), 0, 0],
            chord=ft2m(8, 8),
            twist=0
        ),
        WingXSec(
            airfoil=Airfoil("naca0012"),
            le_loc=[0, 0, ft2m(1)],
            chord=ft2m(3, 8),
            twist=0
        ),
        WingXSec(
            airfoil=Airfoil("naca0012"),
            le_loc=[ft2m(0, 8), 0, ft2m(5)],
            chord=ft2m(2, 8), 
            twist=0
        )
    ]
    vertical_stabilizer = Wing(name="Vertical Stabilizer", xsecs=vs_xsecs,symmetric=false)
    translate!(vertical_stabilizer, [ft2m(16, 11) - ft2m(3, 8), 0, ft2m(-2)])

    # Define the fuselage
    xc =[0,0,ft2m(3),ft2m(5),ft2m(10,4),ft2m(12,4),ft2m(21,11)]
    zc = [ft2m(-1),ft2m(-1),ft2m(-0.85),ft2m(0),ft2m(0.3),ft2m(-0.5,4),ft2m(0.2)]
    radii = [ft2m(0.1), ft2m(1.5), ft2m(1.7), ft2m(2.7), ft2m(2.3),ft2m(1,4), ft2m(0.7)]
    shapes = [2, 3, 7, 7, 7, 5, 3]

    fuse_xsecs = [FuselageXSec(radius=radii[i],xyz_c=[xc[i], 0, zc[i]],shape=shapes[i]) for i in eachindex(xc)]
    fuselage = Fuselage(name="Main Body", xsecs=fuse_xsecs)
    translate!(fuselage, [ft2m(-5), 0, ft2m(-3)])
    # Combine into an airplane
    airplane = Airplane(
        name="Cessna 152",
        wings=[wing, horizontal_stabilizer, vertical_stabilizer],
        fuselages=[fuselage]
    )

    return airplane
end 



 