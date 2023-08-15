abstract type Singularity end

struct LineVortex <: Singularity
    Γ1::Float64
    Γ2::Float64
    x1::Float64
    y1::Float64
    x2::Float64
    y2::Float64
end
