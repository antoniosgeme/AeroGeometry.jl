export cos_space,half_cos_space,rotate2D

cos_space(min,max,npoints) = (max+min)/2 .+ (max-min)/2 * cos.(Array{Float64}(LinRange(π, 0, npoints)))

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
        0         -axis[3]  axis[2];
        axis[3]   0        -axis[1];
       -axis[2]   axis[1]   0
    ]
    
    R = I + sin(θ) * K + (1 - cos(θ)) * (K * K)
    
    # Apply the rotation matrix to the vector
    return R * v
end

cleanup(name) = lowercase(strip(name))
