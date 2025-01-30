function Cessna152()
    ft(feet, inches) = 0.3048 * feet + 0.0254 * inches

    # Define the wings
    wing_xsecs = [
        WingXSec(
            airfoil=Airfoil("naca2412"),
            le_loc=[0, 0, 0],
            chord=ft(5,4),
        ),
        WingXSec(
            airfoil=Airfoil("naca2412"),
            le_loc=[0, ft(7,0), ft(7,0) * sind(1)],
            chord=ft(5,4),
        ),
        WingXSec(
            airfoil=Airfoil("naca0012"),
            le_loc=[ft(4, 3/4) - ft(3, 8 + 1/2), ft(33, 4)/2, ft(33, 4)/2 * sind(1)],
            chord=ft(3, 8 + 1/2),
            twist=0
        )
    ]
    wing = Wing(name="Main Wing", xsecs=wing_xsecs, symmetric=true)

    hs_xsecs = [
        WingXSec(
            airfoil=Airfoil("naca0012"),
            le_loc=[0, 0, 0],
            chord=ft(3,8),
            twist=-2
        ),
        WingXSec(
            airfoil=Airfoil("naca0012"),
            le_loc=[ft(1,0), ft(10,0) / 2, 0],
            chord=ft(2, 4 + 3 / 8),
            twist=-2
        )
    ]
    horizontal_stabilizer = Wing(name="Horizontal Stabilizer", xsecs=hs_xsecs, symmetric=true)
    translate!(horizontal_stabilizer, [4.0648, 0, -0.6096])

    vs_xsecs = [
        WingXSec(
            airfoil=Airfoil("naca0012"),
            le_loc=[ft(-5,0), 0, 0],
            chord=ft(8, 8),
            twist=0
        ),
        WingXSec(
            airfoil=Airfoil("naca0012"),
            le_loc=[ft(0,0), 0, ft(1,0)],
            chord=ft(3, 8),
            twist=0
        ),
        WingXSec(
            airfoil=Airfoil("naca0012"),
            le_loc=[ft(0, 8), 0, ft(5,0)],
            chord=ft(2, 8), 
            twist=0
        )
    ]
    vertical_stabilizer = Wing(name="Vertical Stabilizer", xsecs=vs_xsecs,symmetric=false)
    translate!(vertical_stabilizer, [ft(16, 11) - ft(3, 8), 0, ft(-2,0)])

    # Define the fuselage
    xc =[0,0,ft(3,0),ft(5,0),ft(10,4),ft(12,4),ft(21,11)]
    zc = [ft(-1,0),ft(-1,0),ft(-0.85,0),ft(0,0),ft(0.3,0),ft(-0.5,4),ft(0.2,0)]
    radii = [ft(0.1,0), ft(1.5,0), ft(1.7,0), ft(2.7,0), ft(2.3,0),ft(1,4), ft(0.7,0)]
    shapes = [2, 3, 7, 7, 7, 5, 3]

    fuse_xsecs = [FuselageXSec(radius=radii[i],xyz_c=[xc[i], 0, zc[i]],shape=shapes[i]) for i in eachindex(xc)]
    fuselage = Fuselage(name="Main Body", xsecs=fuse_xsecs)
    translate!(fuselage, [ft(-5,0), 0, ft(-3,0)])
    # Combine into an airplane
    airplane = Airplane(
        name="Cessna 152",
        wings=[wing, horizontal_stabilizer, vertical_stabilizer],
        fuselages=[fuselage]
    )

    return airplane
end 



 