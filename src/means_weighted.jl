"""
    mean(A::AbstractArray, w::AbstractWeights, Estimator::EmpiricalMean)
    
The usual empirical mean estimator. Nothing fancy.
"""
function mean(A::AbstractArray, w::AbstractWeights, Estimator::EmpiricalMean)
    return mean(A, w)
end

"""
    chunk(n, k) 
Divide an set {1,2,...,n} into k blocks of equal size (exept the last one if n is not a multiple of k)
"""
chunk(n, k) = [[(n÷k)*j+1:(n÷k)*(j+1) for j in 0:k-2]; [((n÷k)*(k-1)+1):n]]

"""
MedianOfMean(x::AbstractArray, w::AbstractWeights, k::Integer)

Compute the Median of Mean with `k` groups (it does not permute the samples)
"""
function MoM(x::AbstractArray, w::AbstractWeights, k::Integer)
    keep = findall(!iszero, w)
    x = x[keep]
    w = w[keep]
    n = length(x)
    # @assert n % k == 0 "Check that n=$n is a mutliple of k=$k"
    # m = n ÷ k
    block_mean = zeros(k)
    block_weight = weights(ones(k))
    for (i, block) in enumerate(chunk(n, k))
        block_mean[i] = mean(@view(x[block]), w[block])
        block_weight[i] = sum(w[block]) / sqrt(length(block)) # or sum? or else
    end
    return median(block_mean, block_weight)
end
"""
mean(A::AbstractArray, w::AbstractWeights, Estimator::MoM)

The Median of Mean estimator.
"""
function mean(A::AbstractArray, w::AbstractWeights, Estimator::MedianOfMean)
    return MoM(A, w, Estimator.k)
end

"""
    Z_estimator(x::AbstractArray, w::AbstractWeights, α, ψ::Function; ini = median(x))

Implement Z estimators given a 
"""
function Z_estimator(x::AbstractArray, w::AbstractWeights, α, ψ::Function; ini=median(x))
    #! For some reason it seems that `find_zero` only accepts `Float64`
    f(y) = Rα(y, x, w, α, ψ)
    return find_zero(f, ini)
end

Rα(y, x::AbstractArray, w::AbstractWeights, α, ψ::Function) = sum(w[i] * ψ(α * (x[i] - y)) for i in eachindex(x))

"""
mean(A::AbstractArray, w::AbstractWeights, Estimator::Z_Estimator)

A Z estimator given an influence function `x->ψ(x)` and a scaling parameter `α`.
"""
function mean(A::AbstractArray, w::AbstractWeights, Estimator::Z_Estimator; kwargs...)
    return Z_estimator(A, w, Estimator.α, Estimator.ψ; kwargs...)
end

"""
    LeeVal(x, w, δ; α₀ = 0.0) 

Reference: Optimal Sub-Gaussian Mean Estimation in R by Lee et Valiant 
"""
function LeeVal(x::AbstractArray, w::AbstractWeights, δ; α₀=0.0)
    κ = MoM(x, w, ceil(Int, log(1 / δ)))
    f(α) = sum([min(α * (x[i] - κ)^2, 1) for i in eachindex(x)], w) - 1 / 3 * log(1 / δ)
    α = find_zero(f, α₀)
    return κ + mean([(x[i] - κ) * (1 - min(α * (x[i] - κ)^2, 1)) for i in eachindex(x)], w)
end

"""
    mean(A::AbstractArray, w::AbstractWeights, Estimator::δLeeValiant; kwargs...)

Reference: Optimal Sub-Gaussian Mean Estimation in R by Lee et Valiant 
"""
function mean(A::AbstractArray, w::AbstractWeights, δ::Real, Estimator::LeeValiant; kwargs...)
    return LeeVal(A, w, δ; kwargs...)
end

"""
    MinskerNdaoud(x, w, k, p)

Reference: Robust and efficient mean estimation: an approach based on the properties of self-normalized sums
"""
function MinNda(x, w::AbstractWeights, k, p)
    n = length(x)
    @assert n % k == 0 "Check that n=$n is a mutliple of k=$k"
    m = n ÷ k
    μ̄ = [mean(@view(x[(1+(l-1)*m):(l*m)])) for l in 1:k]
    σ̂p = [stdm(@view(x[(1+(l-1)*m):(l*m)]), μ̄[l], corrected=false) for l in 1:k] .^ p
    return mean(μ̄ ./ σ̂p) * harmmean(σ̂p)
end

"""
    mean(A::AbstractArray, w::AbstractWeights, Estimator::δMinskerNdaoud; kwargs...)

Reference: Robust and efficient mean estimation: an approach based on the properties of self-normalized sums
"""
function mean(A::AbstractArray, w::AbstractWeights, δ::Real, Estimator::MinskerNdaoud; kwargs...)
    k = k_mn(δ)
    return MinNda(A, w, k, Estimator.p; kwargs...)
end