"""
Represents a control surface on an aerodynamic body.

# Fields
- `name::String`: Name of the control surface.
- `xsec_id::Vector{Int64}`: Indices of cross-sections associated with the control surface.
- `deflection::Float64`: Deflection angle of the control surface in degrees (default: 0.0).
- `hinge_point::Float64`: Hinge location along the chord as a fraction (default: 0.75).
- `symmetric::Bool`: Indicates whether the control surface deflects symmetrically (default: true).

# Example
```julia
cs = ControlSurface(name="Aileron", xsec_id=[1, 2, 3], deflection=5.0, hinge_point=0.4, symmetric=false)
```
"""
mutable struct ControlSurface <: AeroComponent
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
        symmetric::Bool = true,
    )
        return new(name, xsec_id, deflection, hinge_point, symmetric)
    end
end

function show(io::IO, cs::ControlSurface)
    println(io, "Control Surface: ", cs.name)
    println(io, "  Deflection: ", cs.deflection, "°")
    println(io, "  Hinge Point: ", cs.hinge_point)
end


"""
Represents a cross-section of a wing, defined by its airfoil shape, leading edge location, chord length, and twist angle.

# Fields
- `airfoil::Airfoil`: The airfoil object representing the shape of this wing section.
- `le_loc::Vector{Float64}`: A 3-element vector specifying the leading edge location `[x, y, z]` of this cross-section in space before any twist is applied.
- `chord::Float64`: The chord length of the wing cross-section.
- `twist::Float64`: The twist angle (in degrees) of this cross-section relative to the root.

# Outer Constructor
The `WingSection` struct has a custom constructor with the following keyword arguments (all optional):
- `airfoil`: The airfoil object (default is `Airfoil("NACA0012")`).
- `le_loc`: The leading edge location before any twist is applied (default is `[0.0, 0.0, 0.0]`).
- `chord`: The chord length (default is `1.0`).
- `twist`: The twist angle (default is `0.0`).

# Example
```julia
xsec = WingSection(le_loc=[1.0, 0.5, 0.0], chord=2.0, twist=5.0)
```
"""
mutable struct WingSection <: AeroComponent
    airfoil::Airfoil
    le_loc::Vector{Float64}
    chord::Float64
    twist::Float64

    # Custom outer constructor with keyword arguments and default values
    function WingSection(;
        airfoil = Airfoil("NACA0012"),
        le_loc = [0.0, 0.0, 0.0],
        chord = 1.0,
        twist = 0.0,
    )
        new(airfoil, le_loc, chord, twist)
    end
end

function show(io::IO, xsec::WingSection)
    println(io, "Wing Cross-Section:")
    println(io, "  Airfoil: ", xsec.airfoil)
    println(io, "  Leading Edge Location: ", xsec.le_loc)
    println(io, "  Chord Length: ", xsec.chord)
    println(io, "  Twist: ", xsec.twist, "°")
end



"""
    Represents an entire wing, consisting of multiple cross-sections and a symmetry property.

    # Fields
    - `name::String`: The name of the wing, useful for identification (e.g., "Main Wing").
    - `sections::Vector{WingSection}`: A vector of `WingSection` objects, each representing a cross-section of the wing. These define the shape and geometry of the wing along its span.
    - `symmetric::Bool`: Indicates whether the wing is symmetric about the y=0 plane:
        - `true`: The wing is mirrored across the centerline (e.g., for a standard airplane wing pair).
        - `false`: The wing is not mirrored (e.g., for a single wing or an asymmetrical design).
    - `control_surfaces::Vector{ControlSurface}`: A vector of `ControlSurface` objects defining movable surfaces on the wing (default: empty vector).

    # Custom Constructor
    The `Wing` struct has a custom constructor with keyword arguments and default values:
    - `name`: The name of the wing (default is `"Unnamed Wing"`).
    - `sections`: A vector of `WingSection` objects (default is an empty vector `Vector{WingSection}()`).
    - `symmetric`: A boolean value indicating symmetry (default is `true`).
    - `control_surfaces`: A vector of `ControlSurface` objects (default is `Vector{ControlSurface}()`).

    # Example
    ```julia
    # Define individual cross-sections of the wing
    xsec1 = WingSection(le_loc=[0.0, 0.0, 0.0], chord=1.0, twist=0.0)
    xsec2 = WingSection(le_loc=[1.0, 0.5, 0.1], chord=0.8, twist=3.0)
    
    # Define a control surface
    cs1 = ControlSurface(name="Aileron", xsec_id=[1, 2], deflection=5.0, hinge_point=0.4, symmetric=false)
    
    # Create the wing with cross-sections and control surfaces
    wing = Wing(name="Example Wing", sections=[xsec1, xsec2], symmetric=true, control_surfaces=[cs1])
    ```
    """
mutable struct Wing <: AeroComponent
    name::String
    sections::Vector{WingSection}
    symmetric::Bool
    control_surfaces::Vector{ControlSurface}

    function Wing(;
        name = "Unnamed Wing",
        sections = Vector{WingSection}(),
        symmetric = true,
        control_surfaces = Vector{ControlSurface}(),
    )
        new(name, sections, symmetric, control_surfaces)
    end
end

function show(io::IO, wing::Wing)
    println(io, "Wing: ", wing.name)
    println(io, "  Number of cross-sections: ", length(wing.sections))
    println(io, "  Symmetric: ", wing.symmetric)
    println(io, "  Number of control surfaces: ", length(wing.control_surfaces))
end


"""
    coordinates(wing::Wing,camberline=false)

Computes the coordinates of each wing cross-section in the global reference frame.
if camberline is true, it returns the wing camberline coordinates instead of the surface coordinates.

"""
function coordinates(wing::Wing; camberline::Bool = false)

    N = length(wing.sections[1].airfoil.x)

    x_surface = zeros(N, length(wing.sections))
    y_surface = zeros(N, length(wing.sections))
    z_surface = zeros(N, length(wing.sections))


    for i = 1:length(wing.sections)
        xsec = wing.sections[i]
        if camberline
            # Get camberline evaluated at the same x-coordinates as the airfoil
            xc = xsec.airfoil.x
            airfoil_coords = camber(xsec.airfoil, xc = xc)
        else
            airfoil_coords = coordinates(xsec.airfoil)
        end
        coords = hcat(
            airfoil_coords[:, 1],
            zeros(size(airfoil_coords[:, 1])),
            airfoil_coords[:, 2],
        )
        chord = xsec.chord
        twist_angle = xsec.twist
        le_loc = xsec.le_loc

        # Compute local wing frame
        xg_local, yg_local, zg_local = compute_frame(wing, i)
        basis = hcat(xg_local, yg_local, zg_local)
        translated_coords = coords * basis .* chord .+ le_loc'

        qc = quarter_chord(xsec.airfoil) .* chord
        quarter_chord_current = le_loc .+ [qc[1], 0, qc[2]]
        twist_angle_rad = deg2rad(twist_angle)

        # Rotate each point around the computed axis
        for (j, point) in enumerate(eachrow(translated_coords))
            rotated_coords =
                rotate_vector(point - quarter_chord_current, yg_local, twist_angle_rad) +
                quarter_chord_current
            x_surface[j, i] = rotated_coords[1]
            y_surface[j, i] = rotated_coords[2]
            z_surface[j, i] = rotated_coords[3]
        end
    end
    return (x_surface, y_surface, z_surface)
end


"""
    get_ghost_idx(wing::Wing)

Figures out where ghost cross sections need to be added to account for control surfaces
"""
function get_ghost_idx(wing::Wing)
    # If control surfaces exist add necessary ghost sections
    ghosts_idx = Vector{Int64}()
    for cs in wing.control_surfaces
        for id in cs.xsec_id
            if id != 1 && id != length(wing.sections) && cs.deflection != 0
                push!(ghosts_idx, id)
            end
        end
    end
    unique!(ghosts_idx)
    return ghosts_idx
end


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
        project_to_YZ(wing.sections[2].le_loc - wing.sections[1].le_loc)
    elseif index == length(wing.sections)
        project_to_YZ(wing.sections[end].le_loc - wing.sections[end-1].le_loc)
    else
        vec_before =
            project_to_YZ(wing.sections[index].le_loc - wing.sections[index-1].le_loc)
        vec_after =
            project_to_YZ(wing.sections[index+1].le_loc - wing.sections[index].le_loc)
        span_vec = normalize(vec_before + vec_after)
        z_scale = sqrt(2 / (1 + dot(vec_before, vec_after)))
        span_vec * z_scale
    end

    zg_local = normalize(cross(xg_local, yg_local))  # Local normal axis

    # Apply twist using a 3D rotation matrix
    twist_angle = deg2rad(wing.sections[index].twist)
    xg_local = rotate_vector(xg_local, yg_local, twist_angle)
    zg_local = rotate_vector(zg_local, yg_local, twist_angle)

    return xg_local, yg_local, zg_local
end

function translate!(wing::Wing, xyz::Vector{Float64})
    for xsec in wing.sections
        xsec.le_loc = xsec.le_loc .+ xyz
    end
end

function deflect_control_surface!(wing::Wing; name = "Aileron", deflection = 0)
    control_surface = get_control_surface(wing, name)
    control_surface.deflection=deflection
    for i in collect(minimum(control_surface.xsec_id):maximum(control_surface.xsec_id))
        deflect_control_surface!(wing.sections[i].airfoil, deflection = deflection)
    end
end


"""
Find a control surface by name in a wing.
"""
function get_control_surface(wing::Wing, name::String)
    cleanup(n) = lowercase(strip(n))
    idx = findfirst(cs -> cleanup(cs.name) == cleanup(name), wing.control_surfaces)
    return isnothing(idx) ? nothing : wing.control_surfaces[idx]
end

function repanel!(wing::Wing, points_per_side)
    for xsec in wing.sections
        repanel!(xsec.airfoil, points_per_side)
    end
    return wing
end


area(xsec::WingSection) = area(xsec.airfoil) * xsec.chord^2


"""
    mesh(wing::Wing; n_span=nothing, n_chord=nothing, camberline=false) -> (points, faces, section_num)

Generate a surface mesh for the given `wing` in the common `(points, faces)` format:

- `points` is an `Npoints × 3` matrix of the 3D vertex coordinates.
- `faces` is an `Nfaces × 4` matrix (for quadrilaterals). Each row contains
  the indices (1-based) into `points` that define one face of the mesh.

If `wing.symmetric == true`, the geometry is mirrored about `y=0`, so the mesh
covers both sides (left and right) of the wing.

# Arguments
- `wing::Wing`: The wing to mesh
- `n_span::Union{Int,Nothing}`: Number of spanwise stations (default: number of wing sections)
- `n_chord::Union{Int,Nothing}`: Number of chordwise points per section (default: number of airfoil points)
- `camberline::Bool`: If true, mesh the camberline instead of the surface

# Returns
- `(points, faces, section_num)`: 
  - `points` is an array of size `(N×M, 3)` where N=n_chord, M=n_span
  - `faces` is an array of size `((N-1)×(M-1), 4)` defining quadrilateral panels
  - `section_num` is an array of size `(N×M,)` indicating the section index for each point
"""
function mesh(
    wing::Wing;
    n_span::Union{Int,Nothing} = nothing,
    n_chord::Union{Int,Nothing} = nothing,
    camberline::Bool = false,
)
    # Get the base coordinates
    x_surf, y_surf, z_surf = coordinates(wing, camberline = camberline)

    N_orig, M_orig = size(x_surf)  # original chordwise and spanwise dimensions

    # Set defaults if not provided
    N = isnothing(n_chord) ? N_orig : n_chord
    M = isnothing(n_span) ? M_orig : n_span

    # If discretization matches original, return the simple mesh
    if N == N_orig && M == M_orig
        return simple_mesh(x_surf, y_surf, z_surf)
    end

    # Otherwise, interpolate to achieve desired discretization
    x_mesh, y_mesh, z_mesh = interpolate_mesh(x_surf, y_surf, z_surf, N, M)

    return simple_mesh(x_mesh, y_mesh, z_mesh)
end

"""
    simple_mesh(x_surf, y_surf, z_surf) -> (points, faces, section_num)

Internal helper to build mesh arrays from coordinate matrices without interpolation.
"""
function simple_mesh(x_surf, y_surf, z_surf)
    N, M = size(x_surf)  # chordwise = N, spanwise = M

    np = N * M
    points = zeros(Float64, np, 3)
    section_num = zeros(Int, np)

    idx = 1
    for j = 1:M
        for i = 1:N
            points[idx, 1] = x_surf[i, j]
            points[idx, 2] = y_surf[i, j]
            points[idx, 3] = z_surf[i, j]
            section_num[idx] = j
            idx += 1
        end
    end

    # Build the faces array for a structured quadrilateral mesh
    nfaces = (N - 1) * (M - 1)
    faces = Array{Int,2}(undef, nfaces, 4)

    face_idx = 1
    for j = 1:(M-1)
        for i = 1:(N-1)
            # Convert (i, j) in [1-based 2D] to the flattened index
            p1 = i + (j - 1) * N
            p2 = i + 1 + (j - 1) * N
            p3 = i + j * N
            p4 = i + 1 + j * N

            # Fill one row of faces with these 4 corner indices
            faces[face_idx, :] .= [p1, p2, p4, p3]
            face_idx += 1
        end
    end

    return points, faces, section_num
end

"""
    interpolate_mesh(x_surf, y_surf, z_surf, N_new, M_new) -> (x_new, y_new, z_new)

Interpolate the surface coordinates to achieve desired chordwise (N_new) and spanwise (M_new) discretization.

This function performs:
1. Spanwise interpolation between wing sections (if M_new != M_orig)
2. Chordwise resampling along each section (if N_new != N_orig)

The interpolation uses linear interpolation in 3D space.
"""
function interpolate_mesh(x_surf, y_surf, z_surf, N_new, M_new)
    N_orig, M_orig = size(x_surf)

    # First, interpolate spanwise if needed
    if M_new != M_orig
        # Create spanwise parameter (0 to 1 along the span)
        η_orig = range(0, 1, length = M_orig)
        η_new = range(0, 1, length = M_new)

        # Interpolate each chordwise station
        x_span = zeros(N_orig, M_new)
        y_span = zeros(N_orig, M_new)
        z_span = zeros(N_orig, M_new)

        for i = 1:N_orig
            # Linear interpolation for each coordinate
            x_span[i, :] = _linear_interp(η_orig, x_surf[i, :], η_new)
            y_span[i, :] = _linear_interp(η_orig, y_surf[i, :], η_new)
            z_span[i, :] = _linear_interp(η_orig, z_surf[i, :], η_new)
        end
    else
        x_span = x_surf
        y_span = y_surf
        z_span = z_surf
    end

    # Then, resample chordwise if needed
    if N_new != N_orig
        x_new = zeros(N_new, M_new)
        y_new = zeros(N_new, M_new)
        z_new = zeros(N_new, M_new)

        # For each spanwise section, resample the chordwise distribution
        for j = 1:M_new
            # Find the min and max x-coordinates (LE and TE)
            x_min = minimum(x_span[:, j])
            x_max = maximum(x_span[:, j])

            # Create uniform distribution in x from LE to TE
            x_new[:, j] = range(x_min, x_max, length = N_new)

            # Interpolate y and z based on x-coordinateå
            y_new[:, j] = _linear_interp(x_span[:, j], y_span[:, j], x_new[:, j])
            z_new[:, j] = _linear_interp(x_span[:, j], z_span[:, j], x_new[:, j])
        end
    else
        x_new = x_span
        y_new = y_span
        z_new = z_span
    end

    return x_new, y_new, z_new
end

"""
    _linear_interp(x, y, x_new) -> y_new

Simple linear interpolation: given points (x, y), evaluate at x_new.
"""
function _linear_interp(x, y, x_new)
    n = length(x_new)
    y_new = zeros(eltype(y), n)

    for i = 1:n
        xi = x_new[i]

        # Find the bracketing indices
        if xi <= x[1]
            y_new[i] = y[1]
        elseif xi >= x[end]
            y_new[i] = y[end]
        else
            # Find j such that x[j] <= xi < x[j+1]
            j = searchsortedlast(x, xi)
            if j >= length(x)
                y_new[i] = y[end]
            else
                # Linear interpolation
                t = (xi - x[j]) / (x[j+1] - x[j])
                y_new[i] = y[j] * (1 - t) + y[j+1] * t
            end
        end
    end

    return y_new
end

"""
    _compute_arc_length(x, y, z) -> s

Compute cumulative arc length along a 3D curve defined by points (x[i], y[i], z[i]).
Returns a vector s where s[i] is the arc length from the start to point i.
"""
function _compute_arc_length(x, y, z)
    n = length(x)
    s = zeros(n)
    s[1] = 0.0

    for i = 2:n
        dx = x[i] - x[i-1]
        dy = y[i] - y[i-1]
        dz = z[i] - z[i-1]
        ds = sqrt(dx^2 + dy^2 + dz^2)
        s[i] = s[i-1] + ds
    end

    return s
end
