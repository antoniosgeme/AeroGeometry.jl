
cos_space(min, max, npoints) =
    (max+min)/2 .+ (max-min)/2 * cos.(Array{Float64}(LinRange(π, 0, npoints)))

half_cos_space(npoints) = 1 .- cos.(LinRange(π, 0, npoints) ./ 2)

function rotate2D(vectors::AbstractMatrix{<:Real}, θ::Real)
    @assert size(vectors, 2) == 2 "Input must be an Nx2 matrix!"
    R = [cosd(θ) -sind(θ); sind(θ) cosd(θ)]
    return vectors * R'
end


"""
    rotate_vector(v, axis, θ)

Rotates a 3D vector `v` around a given axis `axis` by an angle `θ` (in radians).
The axis should be a unit vector.
"""
function rotate_vector(v::Vector{<:Number}, axis::Vector{<:Number}, θ::Number)
    # Ensure the axis is normalized
    axis = normalize(axis)

    # Compute the rotation matrix using the Rodrigues' rotation formula
    K = [
        0 -axis[3] axis[2];
        axis[3] 0 -axis[1];
        -axis[2] axis[1] 0
    ]

    R = I + sin(θ) * K + (1 - cos(θ)) * (K * K)

    # Apply the rotation matrix to the vector
    return R * v
end

cleanup(name) = lowercase(strip(name))


function inside_polygon(xp::AbstractVector, yp::AbstractVector, X, Y)
    n = length(xp)
    @assert n == length(yp) && n ≥ 3
    inside = falses(size(X))

    x1 = xp
    y1 = yp
    x2 = circshift(xp, -1)
    y2 = circshift(yp, -1)

    @inbounds for i = 1:n
        ycross = ((y1[i] .<= Y) .& (Y .< y2[i])) .| ((y2[i] .<= Y) .& (Y .< y1[i]))
        xints = (x2[i] - x1[i]) .* (Y .- y1[i]) ./ (y2[i] - y1[i]) .+ x1[i]
        inside .= xor.(inside, ycross .& (X .< xints))
    end
    return inside
end


function project(vector::Vector{Float64}, plane::Union{Symbol, String})
    plane = Symbol(uppercase(String(plane)))
    if plane == :YZ || plane == :ZY
        projected = [0.0, vector[2], vector[3]]
    elseif plane == :XZ || plane == :ZX
        projected = [vector[1], 0.0, vector[3]]
    elseif plane == :XY || plane == :YX
        projected = [vector[1], vector[2], 0.0]
    else
        error("Invalid plane. Use XY, YX, XZ, ZX, YZ, or ZY")
    end
    return projected 
end