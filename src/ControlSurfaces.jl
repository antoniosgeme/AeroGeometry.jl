"""
Control surface defined as a parametric patch on a Wing.

Coordinates:
- eta ∈ [0,1] root→tip (on one side)
- xi  ∈ [0,1] LE→TE
"""
mutable struct ControlSurface <: AeroComponent
    type::Symbol                  # :Aileron, :Flap, :Elevator, :Rudder, :Slat, etc.

    # Patch bounds
    eta::Tuple{Float64,Float64}   # (eta_in, eta_out)
    xi::Tuple{Float64,Float64}    # (xi_in,  xi_out) 

    # Hinge definition
    xi_hinge::Tuple{Float64,Float64}           # hinge location as fraction of local chord or range

    # Deflection in degrees
    deflection::Float64

    # How it mirrors across y=0 when wing.symmetric == true
    symmetric::Bool

    # Optional scaling for mixers / schedules
    gain::Float64

    function ControlSurface(;
        type=:UnnamedControl,
        eta=(0.0, 1.0),
        xi=(0.75, 1.0),
        xi_hinge=first(xi),
        deflection=0.0,
        symmetric=true,
        gain=1.0,
    )
        # Normalize inputs to tuples for validation
        eta = eta isa Tuple ? eta : (eta, eta)
        xi = xi isa Tuple ? xi : (xi, xi)
        xi_hinge = xi_hinge isa Tuple ? xi_hinge : (xi_hinge, xi_hinge)
        
        # Validate ranges
        0.0 ≤ eta[1] ≤ eta[2] ≤ 1.0 || error("eta must be within [0,1] and ordered.")
        0.0 ≤ xi[1] ≤ xi[2] ≤ 1.0 || error("xi must be within [0,1] and ordered.")
        0.0 ≤ xi_hinge[1] ≤ xi_hinge[2] ≤ 1.0 || error("xi_hinge must be within [0,1].")
        (xi[1] - 1e-12) ≤ xi_hinge[1] && xi_hinge[2] ≤ (xi[2] + 1e-12) || 
            error("xi_hinge range should lie inside xi range.")
        new(type, eta, xi, xi_hinge, deflection, symmetric, gain)
    end
end


function show(io::IO, cs::ControlSurface)
    println(io, "Control Surface: ", cs.type)
    println(io, "  η (root→tip): ", cs.eta)
    println(io, "  ξ range (LE→TE): ", cs.xi)
    println(io, "  Hinge ξ location: ", cs.xi_hinge)
    println(io, "  Deflection (deg): ", cs.deflection)
    println(io, "  Symmetric: ", cs.symmetric)
    println(io, "  Gain: ", cs.gain)
end


function deflect!(cs::ControlSurface, δ::Float64)
    cs.deflection += δ
    return nothing
end