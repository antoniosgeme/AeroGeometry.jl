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