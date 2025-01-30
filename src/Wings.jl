"""
Represents a control surface on an aerodynamic body.

# Fields
- `name::String`: Name of the control surface.
- `xsec_id::Vector{Int64}`: Indices of cross-sections associated with the control surface.
- `deflection::Float64`: Deflection angle of the control surface in degrees (default: 0.0).
- `hinge_point::Float64`: Hinge location along the chord as a fraction (default: 0.5).
- `symmetric::Bool`: Indicates whether the control surface deflects symmetrically (default: true).

# Example
```julia
cs = ControlSurface(name="Aileron", xsec_id=[1, 2, 3], deflection=5.0, hinge_point=0.4, symmetric=false)
```
"""
mutable struct ControlSurface
    name::String  
    xsec_id::Vector{Int64}
    deflection::Float64
    hinge_point::Float64
    symmetric::Bool

    function ControlSurface(; 
        name::String = "Default",
        xsec_id::Vector{Int64} = Int64[],
        deflection::Float64 = 0.0,
        hinge_point::Float64 = 0.75,
        symmetric::Bool = true
    )
        return new(name, xsec_id, deflection, hinge_point, symmetric)
    end
end


"""
Represents a cross-section of a wing, defined by its airfoil shape, leading edge location, chord length, and twist angle.

# Fields
- `airfoil::Airfoil`: The airfoil object representing the shape of this wing section.
- `le_loc::Vector{Float64}`: A 3-element vector specifying the leading edge location `[x, y, z]` of this cross-section in space before any twist is applied.
- `chord::Float64`: The chord length of the wing cross-section.
- `twist::Float64`: The twist angle (in degrees) of this cross-section relative to the root.

# Outer Constructor
The `WingXSec` struct has a custom constructor with the following keyword arguments (all optional):
- `airfoil`: The airfoil object (default is `Airfoil("NACA0012")`).
- `le_loc`: The leading edge location before any twist is applied (default is `[0.0, 0.0, 0.0]`).
- `chord`: The chord length (default is `1.0`).
- `twist`: The twist angle (default is `0.0`).

# Example
```julia
xsec = WingXSec(le_loc=[1.0, 0.5, 0.0], chord=2.0, twist=5.0)
```
"""
mutable struct WingXSec
    airfoil::Airfoil
    le_loc::Vector{Float64}
    chord::Float64
    twist::Float64

    # Custom outer constructor with keyword arguments and default values
    function WingXSec(; airfoil=Airfoil("NACA0012"), le_loc=[0.0, 0.0, 0.0], chord=1.0, twist=0.0)
        new(airfoil, le_loc, chord, twist)
    end
end 


"""
    Represents an entire wing, consisting of multiple cross-sections and a symmetry property.

    # Fields
    - `name::String`: The name of the wing, useful for identification (e.g., "Main Wing").
    - `xsecs::Vector{WingXSec}`: A vector of `WingXSec` objects, each representing a cross-section of the wing. These define the shape and geometry of the wing along its span.
    - `symmetric::Bool`: Indicates whether the wing is symmetric about the y=0 plane:
        - `true`: The wing is mirrored across the centerline (e.g., for a standard airplane wing pair).
        - `false`: The wing is not mirrored (e.g., for a single wing or an asymmetrical design).
    - `control_surfaces::Vector{ControlSurface}`: A vector of `ControlSurface` objects defining movable surfaces on the wing (default: empty vector).

    # Custom Constructor
    The `Wing` struct has a custom constructor with keyword arguments and default values:
    - `name`: The name of the wing (default is `"Unnamed Wing"`).
    - `xsecs`: A vector of `WingXSec` objects (default is an empty vector `Vector{WingXSec}()`).
    - `symmetric`: A boolean value indicating symmetry (default is `true`).
    - `control_surfaces`: A vector of `ControlSurface` objects (default is `Vector{ControlSurface}()`).

    # Example
    ```julia
    # Define individual cross-sections of the wing
    xsec1 = WingXSec(le_loc=[0.0, 0.0, 0.0], chord=1.0, twist=0.0)
    xsec2 = WingXSec(le_loc=[1.0, 0.5, 0.1], chord=0.8, twist=3.0)
    
    # Define a control surface
    cs1 = ControlSurface(name="Aileron", xsec_id=[1, 2], deflection=5.0, hinge_point=0.4, symmetric=false)
    
    # Create the wing with cross-sections and control surfaces
    wing = Wing(name="Example Wing", xsecs=[xsec1, xsec2], symmetric=true, control_surfaces=[cs1])
    ```
    """
mutable struct Wing
    name::String
    xsecs::Vector{WingXSec}
    symmetric::Bool
    control_surfaces::Vector{ControlSurface}

    function Wing(; name="Unnamed Wing", xsecs=Vector{WingXSec}(), symmetric=true, control_surfaces=Vector{ControlSurface}())
        new(name, xsecs, symmetric, control_surfaces)
    end
end

"""
    coordinates(wing::Wing)

Computes the coordinates of each wing cross-section in the global reference frame
"""
function coordinates(wing::Wing)
    N = length(wing.xsecs[1].airfoil.coordinates[:,1])
    x_surface = zeros(N,length(wing.xsecs))
    y_surface = zeros(N,length(wing.xsecs))
    z_surface = zeros(N,length(wing.xsecs))

    for i in 1:length(wing.xsecs)
        xsec = wing.xsecs[i]
        coords = hcat(xsec.airfoil.coordinates[:,1], zeros(size(xsec.airfoil.coordinates[:,1])), xsec.airfoil.coordinates[:,2])
        chord = xsec.chord
        twist_angle = xsec.twist
        le_loc = xsec.le_loc

        # Compute local wing frame
        xg_local, yg_local, zg_local = compute_frame(wing, i)
        basis = hcat(xg_local, yg_local, zg_local)
        translated_coords = coords * basis .* chord .+ le_loc'

        

       
        qc = quarter_chord(xsec.airfoil) .* chord
        quarter_chord_current = le_loc .+ [qc[1],0,qc[2]]
        twist_angle_rad = deg2rad(twist_angle)

        # Rotate each point around the computed axis
        for (j, point) in enumerate(eachrow(translated_coords))
            rotated_coords = rotate_vector(point - quarter_chord_current, yg_local, twist_angle_rad) + quarter_chord_current
            x_surface[j,i] = rotated_coords[1]
            y_surface[j,i] = rotated_coords[2]
            z_surface[j,i] = rotated_coords[3]
        end 
    end
    return (x_surface, y_surface, z_surface)
end



area(xsec::WingXSec) = area(xsec.airfoil) * xsec.chord^2



function compute_frame(wing::Wing, index::Int)
    """
    Computes the local reference frame for a specific cross-section (XSec) of the wing.

    Args:
        wing::Wing: Wing object containing cross-sections.
        index::Int: Index of the cross-section.

    Returns:
        Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}: (xg_local, yg_local, zg_local),
        the local reference frame axes as unit vectors (chordwise, spanwise, normal).
    """

    # Helper: Project a vector onto the YZ plane and normalize it
    project_to_YZ(vector::Vector{Float64}) = [0.0, vector[2], vector[3]] / norm(vector[2:3])

    xg_local = [1.0, 0.0, 0.0]  # Local chordwise axis

    # Determine the spanwise (yg_local) direction
    yg_local = if index == 1
        project_to_YZ(wing.xsecs[2].le_loc - wing.xsecs[1].le_loc)
    elseif index == length(wing.xsecs)
        project_to_YZ(wing.xsecs[end].le_loc - wing.xsecs[end-1].le_loc)
    else
        vec_before = project_to_YZ(wing.xsecs[index].le_loc - wing.xsecs[index-1].le_loc)
        vec_after = project_to_YZ(wing.xsecs[index+1].le_loc - wing.xsecs[index].le_loc)
        span_vec = normalize(vec_before + vec_after)
        z_scale = sqrt(2 / (1 + dot(vec_before, vec_after)))
        span_vec * z_scale
    end

    zg_local = normalize(cross(xg_local, yg_local))  # Local normal axis

    # Apply twist using a 3D rotation matrix
    twist_angle = deg2rad(wing.xsecs[index].twist)
    xg_local = rotate_vector(xg_local, yg_local, twist_angle) 
    zg_local = rotate_vector(zg_local, yg_local, twist_angle) 

    return xg_local, yg_local, zg_local
end

function translate!(wing::Wing, xyz::Vector{Float64})
    for xsec in wing.xsecs
        xsec.le_loc = xsec.le_loc .+ xyz
    end 
end




