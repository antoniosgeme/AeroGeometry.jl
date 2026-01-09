
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
