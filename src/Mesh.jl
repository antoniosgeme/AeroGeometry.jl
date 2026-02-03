function mesh(
    wing::Wing;
    spanwise::Int = 20,
    chordwise::Int = 10,
    spanwise_spacing::Union{String, Symbol} = "uniform",
    chordwise_spacing::Union{String, Symbol} = "uniform",
    camberline::Bool = true,
)
    X, Y, Z = coordinates(wing, camberline = camberline)
    N_orig, M_orig = size(X)

    # Spanwise parameter
    η_orig = range(0.0, 1.0, length = M_orig)
    spacing = lowercase_symbol(spanwise_spacing)
    η_new = if spacing === :cosine || spacing === :cos
        0.5 .- 0.5 .* cos.(range(0, π, length = spanwise))
    elseif spacing === :uniform
        range(0.0, 1.0, length = spanwise)
    else
        error("Invalid spanwise spacing option. Use 'uniform' or 'cosine'.")
    end

    # Chordwise parameter (0..1)
    ξ_orig = range(0.0, 1.0, length = N_orig)
    spacing = lowercase_symbol(chordwise_spacing)
    ξ_new = if spacing === :cosine || spacing === :cos
        0.5 .- 0.5 .* cos.(range(0, π, length = chordwise))
    elseif spacing === :uniform
        range(0.0, 1.0, length = chordwise)
    else
        error("Invalid chordwise spacing option. Use 'uniform' or 'cosine'.")
    end

    k_span = min(3, M_orig - 1)
    k_chord = min(3, N_orig - 1)

    # Spanwise interpolation for each chordwise index
    Xs = zeros(N_orig, spanwise)
    Ys = zeros(N_orig, spanwise)
    Zs = zeros(N_orig, spanwise)
    for i = 1:N_orig
        for (k, η) in enumerate(η_new)
            if η <= η_orig[1]
                Xs[i, k] = X[i, 1]; Ys[i, k] = Y[i, 1]; Zs[i, k] = Z[i, 1]
            elseif η >= η_orig[end]
                Xs[i, k] = X[i, end]; Ys[i, k] = Y[i, end]; Zs[i, k] = Z[i, end]
            else
                j = searchsortedlast(η_orig, η)
                t = (η - η_orig[j]) / (η_orig[j + 1] - η_orig[j])
                Xs[i, k] = X[i, j] * (1 - t) + X[i, j + 1] * t
                Ys[i, k] = Y[i, j] * (1 - t) + Y[i, j + 1] * t
                Zs[i, k] = Z[i, j] * (1 - t) + Z[i, j + 1] * t
            end
        end
    end

     # Chordwise resampling
    Xc = zeros(chordwise, spanwise)
    Yc = zeros(chordwise, spanwise)
    Zc = zeros(chordwise, spanwise)

    for j = 1:spanwise
        if camberline
            # Arc-length parameter for camberline
            s = zeros(N_orig)
            for i = 2:N_orig
                dx = Xs[i, j] - Xs[i - 1, j]
                dy = Ys[i, j] - Ys[i - 1, j]
                dz = Zs[i, j] - Zs[i - 1, j]
                s[i] = s[i - 1] + sqrt(dx^2 + dy^2 + dz^2)
            end
            t_new = s[end] .* ξ_new

            sx = Spline1D(s, Xs[:, j], k = k_chord)
            sy = Spline1D(s, Ys[:, j], k = k_chord)
            sz = Spline1D(s, Zs[:, j], k = k_chord)

            Xc[:, j] = evaluate(sx, t_new)
            Yc[:, j] = evaluate(sy, t_new)
            Zc[:, j] = evaluate(sz, t_new)
        else
            # Index-based parameter for closed surface
            sx = Spline1D(ξ_orig, Xs[:, j], k = k_chord)
            sy = Spline1D(ξ_orig, Ys[:, j], k = k_chord)
            sz = Spline1D(ξ_orig, Zs[:, j], k = k_chord)

            Xc[:, j] = evaluate(sx, ξ_new)
            Yc[:, j] = evaluate(sy, ξ_new)
            Zc[:, j] = evaluate(sz, ξ_new)
        end
    end

    return Xc, Yc, Zc
end


function convert_to_unstructured(
    X::Matrix{Float64},
    Y::Matrix{Float64},
    Z::Matrix{Float64},
)
    N, M = size(X)
    points = hcat(vec(X), vec(Y), vec(Z))

    # Create connectivity (triangulated quads)
    n_quads = (N - 1) * (M - 1)
    faces = Matrix{Int}(undef, 2 * n_quads, 3)
    idx = 1
    for j in 1:(M - 1)
        for i in 1:(N - 1)
            n1 = i + (j - 1) * N
            n2 = n1 + 1
            n3 = n2 + N
            n4 = n1 + N
            faces[idx, :] .= (n1, n2, n3)
            faces[idx + 1, :] .= (n1, n3, n4)
            idx += 2
        end
    end

    return points, faces
end