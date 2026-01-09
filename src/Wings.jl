"""
Represents a cross-section of a wing, defined by its airfoil, leading edge location in the wing-fixed axis, chord length, and twist angle.

# Fields
- `airfoil::Airfoil`: The airfoil object representing the shape of this wing section.
- `position::Vector{Float64}`: A 3-element vector specifying the leading edge location `[x, y, z]` of this cross-section in the wing coordinate system before any twist is applied.
- `chord::Float64`: The chord length of the wing cross-section.
- `twist::Float64`: The twist angle (in degrees) of this cross-section relative to the root.

# Example
```julia
xsec = WingSection(position=[1.0, 0.5, 0.0], chord=2.0, twist=5.0)
```

# Notes
- The leading edge location is defined in the wing coordinate system and applying the twist. 
"""
mutable struct WingSection <: AeroComponent
    airfoil::Airfoil
    position::Vector{Float64}
    chord::Float64
    twist::Float64

    # Custom outer constructor with keyword arguments and default values
    function WingSection(;
        airfoil = Airfoil("NACA0012"),
        position = [0.0, 0.0, 0.0],
        chord = 1.0,
        twist = 0.0,
    )
        new(airfoil, position, chord, twist)
    end
end

function show(io::IO, xsec::WingSection)
    println(io, "Wing Cross-Section:")
    println(io, "  Airfoil: ", xsec.airfoil)
    println(io, "  Leading Edge Location: ", xsec.position)
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

    # Notes
    - The cross-sections in `sections` are automatically sorted by their y-coordinate of the leading edge location upon creation of the `Wing` object.

    # Example
    ```julia
    # Define individual cross-sections of the wing
    xsec1 = WingSection(position=[0.0, 0.0, 0.0], chord=1.0, twist=0.0)
    xsec2 = WingSection(position=[1.0, 0.5, 0.1], chord=0.8, twist=3.0)
    
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
    position::Vector{Float64}  
    orientation::Matrix{Float64} 

    function Wing(;
        name = "Unnamed Wing",
        sections = [WingSection(), WingSection(position=[0, 1.,0])],
        symmetric = true,
        control_surfaces = Vector{ControlSurface}(),
        position = [0.0, 0.0, 0.0],
        orientation = Matrix{Float64}(I, 3, 3)

    )
        
        # reorder section by y-coordinate of leading edge
        sections = sort(sections, by = xsec -> xsec.position[2])
        # warn if wing is symmetric but has at least one section with positive and negative y-coordinate
        if symmetric
            has_positive = any(xsec -> xsec.position[2] > 0, sections)
            has_negative = any(xsec -> xsec.position[2] < 0, sections)
            if has_positive && has_negative
                @warn "Wing is symmetric but has sections with both positive and negative y-coordinates. This may lead to incorrect geometry."
            end
        end

        new(name, sections, symmetric, control_surfaces, position, orientation)
    end
end

function show(io::IO, wing::Wing)
    println(io, "Wing: ", wing.name)
    println(io, "  Number of cross-sections: ", length(wing.sections))
    println(io, "  Symmetric: ", wing.symmetric)
    println(io, "  Number of control surfaces: ", length(wing.control_surfaces))
end


"""
    hub_index(wing::Wing)

Returns the index of the hub section (at y=0) for symmetric wings, or 0 if wing is not symmetric or no hub section exists.
"""
function hub_index(wing::Wing)
    idx = wing.symmetric ? something(findfirst(xsec -> xsec.position[2] == 0, wing.sections), 0) : 0
    if idx > 1 && wing.symmetric
        @warn "Hub section found at index $idx which is not the first section. This may lead to unexpected geometry."
    end
    return idx
end    

"""
    coordinates(wing::Wing,camberline=false)

Computes the coordinates of each wing cross-section in the global coordinate system.
if camberline is true, it returns the wing camberline coordinates instead of the surface coordinates.

# Notes 
- If a wing is symmetric and its hub section is at y=0, that section is projected onto the y=0 plane. to avoid self-intersections.
"""
function coordinates(wing::Wing; camberline::Bool = false)
    
    # find max N of all sections 
    N = maximum(length(xsec.airfoil.x) for xsec in wing.sections)
    repanel!(wing, N)
    M = length(wing.sections)

    X = zeros(N, M)
    Y = zeros(N, M)
    Z = zeros(N, M)

    hub = hub_index(wing)
    
    for (i, xsec) in enumerate(wing.sections)
        if camberline
            xc = range(0, 1, length=N)
            airfoil_coords = camber(xsec.airfoil, xc=xc)
        else
            airfoil_coords = coordinates(xsec.airfoil)
        end
        coords = hcat(
            airfoil_coords[:, 1],
            zeros(size(airfoil_coords[:, 1])),
            airfoil_coords[:, 2],
        )
        chord = xsec.chord
        twist_rad = deg2rad(xsec.twist)
        position = xsec.position


        xg_local, yg_local, zg_local = compute_frame(wing, i)
        basis = hcat(xg_local, yg_local, zg_local)
        if hub == i
            # Project onto XY plane by zeroing the Z component and renormalizing
            basis[:, 1] = normalize([basis[1, 1], basis[2, 1], 0.0])  # xg_local projected
            basis[:, 2] = normalize([basis[1, 2], basis[2, 2], 0.0])  # yg_local projected
            basis[:, 3] = [0.0, 0.0, 1.0]  # zg_local becomes pure Z
        end
        translated_coords = coords * basis' .* chord .+ position'
        qc = quarter_chord(xsec) 

        # Rotate each point around the computed axis
        for (j, point) in enumerate(eachrow(translated_coords))
            rotated_coords = rotate_vector(point - qc, yg_local, twist_rad) + qc
            X[j, i] = rotated_coords[1]
            Y[j, i] = rotated_coords[2]
            Z[j, i] = rotated_coords[3]
        end
    end

    if wing.symmetric
        mirror_idx = filter(i -> i != hub, M:-1:1)
        X = hcat(X[:, mirror_idx], X)
        Y = hcat(-Y[:, mirror_idx], Y)
        Z = hcat(Z[:, mirror_idx], Z)
    end

    X_global = similar(X)
    Y_global = similar(Y)
    Z_global = similar(Z)
    for i in 1:size(X, 1), j in 1:size(X, 2)
        local_point = [X[i,j], Y[i,j], Z[i,j]]
        global_point = wing.orientation * local_point + wing.position
        X_global[i,j] = global_point[1]
        Y_global[i,j] = global_point[2]
        Z_global[i,j] = global_point[3]
    end

    return (X_global, Y_global, Z_global)
end

coordinates(section::WingSection) = coordinates(section.airfoil)


"""
Computes the local reference frame for a specific cross-section of the wing.

Args:
    wing::Wing: Wing object containing cross-sections.
    index::Int: Index of the cross-section.

Returns:
    Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}: (xg_local, yg_local, zg_local),
    the local reference frame axes as unit vectors (chordwise, spanwise, normal).
"""
function compute_frame(wing::Wing, index::Int)
    # Helper: Project a vector onto the YZ plane and normalize it
    xg_local = [1.0, 0.0, 0.0]  # Local chordwise axis

    # Determine the spanwise (yg_local) direction
    yg_local = if index == 1
        normalize(project(wing.sections[2].position - wing.sections[1].position, :YZ))
    elseif index == length(wing.sections)
        normalize(project(wing.sections[end].position - wing.sections[end-1].position, :YZ))
    else
        vec_before = normalize(project(wing.sections[index].position - wing.sections[index-1].position, :YZ))
        vec_after = normalize(project(wing.sections[index+1].position - wing.sections[index].position, :YZ))
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

"""
    translate!(wing::Wing, xyz::Vector{Float64})

Translate the entire wing in global coordinates.
"""
function translate!(wing::Wing, xyz::Vector{<:Real})
    wing.position .+= xyz
    return wing
end


"""
    rotate!(wing::Wing; axis::Vector{<:Real}=[1,0,0], angle::Real=0.0, pivot::Vector{<:Real}=[0,0,0])

Rotate the entire wing around a specified axis by a given angle.

# Arguments
- `wing::Wing`: The wing to rotate (modified in-place)
- `axis::Vector{<:Real}=[1,0,0]`: The rotation axis as a 3D vector in global coordinates. Default is the X-axis.
- `angle::Real=0.0`: The rotation angle in degrees
- `pivot::Vector{<:Real}=[0,0,0]`: The pivot point for rotation in local wing coordinates. Default is the origin.

# Notes
- The pivot is expressed in local coordinates to allow for intuitive rotations around points on the wing.

# Examples
```julia
# Rotate wing 45° around Y-axis at origin
rotate!(wing, axis=[0,1,0], angle=45.0)

# Rotate wing 30° around Z-axis at a specific point
rotate!(wing, axis=[0,0,1], angle=30.0, pivot=[0, 3, 0])
```
"""
function rotate!(wing::Wing; axis::Vector{<:Real}=[1,0,0], angle::Real=0.0, pivot::Vector{<:Real}=[0,0,0])
    pivot_global = wing.orientation * pivot + wing.position
    axis_normalized = normalize(axis)
    # Rodrigues rotation formula
    K = [0 -axis_normalized[3] axis_normalized[2];
         axis_normalized[3] 0 -axis_normalized[1];
         -axis_normalized[2] axis_normalized[1] 0]
    R = I + sind(angle) * K + (1 - cosd(angle)) * K^2
    wing.position = R * (wing.position - pivot_global) + pivot_global
    wing.orientation = R * wing.orientation
    return wing
end


function repanel!(wing::Wing, n::Int)
    for xsec in wing.sections
        repanel!(xsec.airfoil, n)
    end
    return wing
end

function repanel(wing::Wing, n::Int)
    wing_copy = deepcopy(wing)
    repanel!(wing_copy, n)
    return wing_copy
end

"""
    quarter_chord(section::WingSection)

Get the quarter chord point of a wing section in 3D space.
"""
function quarter_chord(section::WingSection) 
    # Get quarter chord of 3D wing section 
    qc = quarter_chord(section.airfoil) .* section.chord
    qc3d = section.position .+ [qc[1], 0, qc[2]]
    return qc3d
end 



"""
    span(wing::Wing; centerline::Bool = true)

Compute the span of a wing by measuring the distance between quarter chord points
of wing sections.

# Arguments
- `wing::Wing`: The wing to measure
- `centerline::Bool=true`: If true and wing is symmetric, includes the distance across 
  the two innermost wing sections

# Returns
- `Float64`: The span in the same units as the wing geometry
"""
function span(wing::Wing; centerline::Bool = false)
    qc_points = [quarter_chord(xsec) for xsec in wing.sections]
    
    semi_span = 0.0
    for i in 1:(length(qc_points)-1)
        semi_span += norm(qc_points[i+1] - qc_points[i])
    end

    if wing.symmetric
        full_span = 2*semi_span
        if centerline 
            d = abs(qc_points[1][2])
            full_span += 2d
        end 
    else
        full_span = semi_span
    end

    return full_span
end

area(xsec::WingSection) = area(xsec.airfoil) * xsec.chord^2


"""
    area(wing::Wing; type::Symbol=:planform, centerline::Bool=false)

Computes the wing area.

For symmetric wings, both left/right sides are included in the total area.

# Arguments
- `wing::Wing`: The wing to measure
- `type`: Either `:planform` (mean camber surface area) or `:wetted` (actual surface area including top/bottom). Can be provided as a symbol or string.
- `centerline::Bool=false`: If true and wing is symmetric with root offset from y=0,
    includes fictitious rectangular area from innermost section to centerline
"""
function area(wing::Wing; type::Union{Symbol, String}=:planform, centerline::Bool=false)
    type_symbol = Symbol(lowercase(String(type)))
    if type_symbol == :planform
        X,Y,Z = coordinates(wing, camberline=true)
    elseif type_symbol == :wetted
        X,Y,Z = coordinates(wing, camberline=false)
    end 
    
    y_offset, root_idx = findmin(xsec -> abs(quarter_chord(xsec)[2]), wing.sections)
    if y_offset > 1e-10 && !centerline && wing.symmetric
       center_idx = root_idx + 1
    else
        center_idx = 0
    end
        
    A = 0.0
    N, M = size(X)
    for j in 1:M-1
        if j !== center_idx
            for i in 1:N-1
                # Define the four corners of the quad panel
                r00 = [X[i, j], Y[i, j], Z[i, j]]
                r10 = [X[i+1, j], Y[i+1, j], Z[i+1, j]]
                r01 = [X[i, j+1], Y[i, j+1], Z[i, j+1]]
                r11 = [X[i+1, j+1], Y[i+1, j+1], Z[i+1, j+1]]

                # Triangle 1: (r00, r10, r01)
                a = r10 - r00
                b = r01 - r00
                A += 0.5 * norm(cross(a, b))

                # Triangle 2: (r10, r11, r01)
                a = r11 - r10
                b = r01 - r10
                A += 0.5 * norm(cross(a, b))
            end 
        end 
    end
    return A
end

"""
    area(sections::Vector{WingSection})
Computes the areas of individual wing sections defined by consecutive WingSection objects.
# Arguments
- `sections::Vector{WingSection}`: A vector of WingSection objects defining the wing sections.
- `type::Symbol`: Either `:planform` (mean camber surface area) or `:wetted` (actual surface area including top/bottom). Can be provided as a symbol or string.
# Returns
- `Vector{Float64}`: A vector containing the area of each wing section.
"""
function area(sections::Vector{WingSection} ;type::Union{Symbol, String}=:planform)
    areas = zeros(Float64, length(sections)-1)
    for i in 1:(length(sections)-1)
        inner = sections[i]
        outer = sections[i+1]
        temp_wing = Wing(sections=[inner, outer])
        section_area = area(temp_wing, type=type)
        areas[i] = section_area
    end
    return areas
end

aspect_ratio(wing::Wing; centerline::Bool=false) = span(wing; centerline=centerline)^2 / area(wing; centerline=centerline)

mean_geometric_chord(wing::Wing) = area(wing) / span(wing)

function mean_aerodynamic_chord(wing::Wing)
    MAC_lengths = zeros(Float64, length(wing.sections)-1)
    areas = area(wing.sections)
    
    for i in 1:(length(wing.sections)-1)
        c1 = wing.sections[i].chord
        c2 = wing.sections[i+1].chord
        λ = c2 / c1
        d = (2/3) * c1 * ((1 + λ + λ^2) / (1 + λ))
        MAC_lengths[i] = d
    end
    MAC = sum(MAC_lengths .* areas) / sum(areas)
    return MAC
end

    

"""
    aerodynamic_center(wing::Wing; chord_fraction::Float64=0.25)

Computes the location of the aerodynamic center of the wing.

# Arguments
- `wing::Wing`: The wing object
- `chord_fraction::Float64`: The position of the aerodynamic center along the MAC, as a fraction of MAC length (default: 0.25)

# Notes
- Typically, `chord_fraction` is 0.25 for a subsonic wing
"""
function aerodynamic_center(wing::Wing; chord_fraction::Float64=0.25)
    aerocenters = Vector{Vector{Float64}}(undef, length(wing.sections)-1)
    areas = zeros(Float64, length(wing.sections)-1)
    
    for i in 1:(length(wing.sections)-1)
        inner = wing.sections[i]
        outer = wing.sections[i+1]
        c1 = inner.chord
        c2 = outer.chord
        d1 = inner.position
        d2 = outer.position
        
        λ = c2 / c1
        MAC = (2/3) * c1 * (1 + λ + λ^2) / (1 + λ)
        MAC_le = d1 + (d2 - d1) * (1 + 2 * λ) / (3 + 3 * λ)
        
        ψ = deg2rad((inner.twist + outer.twist) / 2)
        
        chord_vector = [chord_fraction * MAC, 0.0, 0.0]
        yg_local = normalize(project(d2 - d1, :YZ))
        rotated_chord = rotate_vector(chord_vector, yg_local, ψ)
        AC = MAC_le + rotated_chord
        
        # Calculate section area using temporary wing
        temp_wing = Wing(sections=[inner, outer])
        section_area = area(temp_wing)
        
        aerocenters[i] = AC
        areas[i] = section_area
    end
    
    ac = sum(aerocenters .* areas) / sum(areas)
    
    # If wing is symmetric, force y-coordinate to zero
    if wing.symmetric
        ac[2] = 0.0
    end
    
    return ac
end
