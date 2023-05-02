module tools

export cos_space,half_cos_space,rotate2D

cos_space(min,max,npoints) = (max+min)/2 .+ (max-min)/2 * cos.(Array{Float64}(LinRange(π, 0, npoints)))

half_cos_space(npoints) = 1 .- cos.(LinRange(π, 0, npoints) ./ 2)

rotate2D(v,θ) = [ cos(θ) -sin(θ)
                  sin(θ)  cos(θ) ] * v

end
