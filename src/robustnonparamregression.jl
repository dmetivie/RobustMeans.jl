function r̂ₙ(x, X, Y, kernel)
    Nₕ = [kernel.f((X[i] - x) / kernel.h) for i in eachindex(X)]
    sumNₕ = sum(Nₕ)
    if sumNₕ == 0
        return 0
    else
        return sum(Nₕ[i] * Y[i] for i in eachindex(Y)) / sum(Nₕ)
    end
end

function r̂ₙMOM(x, X, Y, kernel, m)
    N = length(X)
    mom = zeros(m)
    for (i, c) in enumerate(RobustMeans.chunk(N, m))
        mom[i] = r̂ₙ(x, X[c], Y[c], kernel)
    end
    return median(mom)
end