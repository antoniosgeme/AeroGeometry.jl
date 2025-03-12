mutable struct Airplane
    name::String
    fuselages::Vector{Fuselage}
    wings::Vector{Wing}

    function Airplane(;
        name::String = "Untitled",
        fuselages::Vector{Fuselage} = Fuselage[],
        wings::Vector{Wing} = Wing[]
    )
        new(name, fuselages, wings)
    end
end


function show(io::IO, airplane::Airplane)
    println(io, "Airplane: ", airplane.name)
    println(io, "  Number of fuselages: ", length(airplane.fuselages))
    println(io, "  Number of wings: ", length(airplane.wings))
end

function deflect_control_surface!(airplane::Airplane, deflections::Dict{String, <:Number})

    # Clean the deflection dictionary keys
    cleaned_deflections = Dict(cleanup(k) => v for (k, v) in deflections)

    for wing in airplane.wings
        for cs in wing.control_surfaces
            cs_name = cleanup(cs.name)
            if haskey(cleaned_deflections, cs_name)
                deflect_control_surface!(wing, name=cs.name, deflection=cleaned_deflections[cs_name])
            end
        end
    end
end


"""
Find a control surface by name in an entire airplane.
"""
function get_control_surface(airplane::Airplane, name::String)
    cleanup(n) = lowercase(strip(n))
    for wing in airplane.wings
        cs = get_control_surface(wing, name)
        if cs !== nothing
            return cs
        end
    end
    return nothing
end