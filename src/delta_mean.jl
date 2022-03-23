abstract type δ_MeanEstimator <: MeanEstimator end

k_mom(δ) = ceil(Int, 8log(1 / δ))
struct δMoM{T} <: δ_MeanEstimator
    δ::T
end

"""
    mean(A::AbstractArray, MeanEstimator::δMoM)

"""
function mean(A::AbstractArray, MeanEstimator::δMoM)
    k = k_mom(MeanEstimator.δ)
    return mean(A, MoM(k))
end

"""
    TrimmedMean(z, ε) 

"""
function TrimmedMean(z, ε)
    n2 = length(z)
    @assert ε < 1 "ε > 1: you cannot remove more"
    @assert iseven(n2)
    n = n2 ÷ 2
    y = z[1:n]
    x = z[n+1:n2]
    sort!(y)
    i_min = ceil(Int, ε * n)
    i_max = n - floor(Int, ε * n)
    α = y[i_min]
    β = y[i_max]
    return TrimmedMean(x, α, β)
end
ε_trimmed(δ, n) = 32log(8 / δ) / (3n)

"""
     mean(A::AbstractArray, MeanEstimator::δTrimmedMean)

"""
struct δTrimM{T} <: δ_MeanEstimator
    δ::T
end
function mean(A::AbstractArray, MeanEstimator::δTrimM)
    ε = ε_trimmed(MeanEstimator.δ, length(A))
    return TrimmedMean(A, ε)
end

struct δCatoni{T,F} <: δ_MeanEstimator
    δ::T
    σ::F
end
α_Catoni(δ, n, σ) = sqrt((2log(2 / δ)) / (n * (1 + (2log(2 / δ)) / (n - 2log(2 / δ))))) / σ
"""
   mean(A::AbstractArray, MeanEstimator::δCatoni, kwargs...)

Reference: Catoni 
"""
function mean(A::AbstractArray, MeanEstimator::δCatoni, kwargs...)
    n = length(A)
    α = α_Catoni(MeanEstimator.δ, n, MeanEstimator.σ)
    z = Z_Estimator(α, ψ_Catoni)
    return mean(A, z; kwargs...)
end

struct δHuber{T,F} <: δ_MeanEstimator
    δ::T
    σ::F
end
α_Huber(δ, n, σ) = sqrt(2log(4 / δ) / n) / σ #sqrt(2 * (log(2 / δ) - log(1 - exp(-n / 2) / δ)) / n) / σ
"""
   mean(A::AbstractArray, MeanEstimator::δCatoni, kwargs...)

Reference: Huber 
"""
function mean(A::AbstractArray, MeanEstimator::δHuber, kwargs...)
    n = length(A)
    α = α_Huber(MeanEstimator.δ, n, MeanEstimator.σ)
    z = Z_Estimator(α, ψ_Huber)
    mean(A, z; kwargs...)
end

"""
    LeeValiant(x, δ; α₀ = 0.0) 

Reference: Optimal Sub-Gaussian Mean Estimation in R by Lee et Valiant 
"""
function LeeValiant(x, δ; α₀ = 0.0)
    κ = MedianOfMean(x, ceil(Int, log(1 / δ)))
    f(α) = sum(min(α * (x[i] - κ)^2, 1) for i in eachindex(x)) - 1 / 3 * log(1 / δ)
    α = find_zero(f, α₀)
    return κ + mean((x[i] - κ) * (1 - min(α * (x[i] - κ)^2, 1)) for i in eachindex(x))
end
struct δLeeValiant{T} <: δ_MeanEstimator
    δ::T
end

"""
    mean(A::AbstractArray, MeanEstimator::δLeeValiant; kwargs...)

Reference: Optimal Sub-Gaussian Mean Estimation in R by Lee et Valiant 
"""
function mean(A::AbstractArray, MeanEstimator::δLeeValiant; kwargs...)
    return LeeValiant(A, MeanEstimator.δ; kwargs...)
end

"""
    MinskerNdaoud(x, k, p)

Reference: Robust and efficient mean estimation: an approach based on the properties of self-normalized sums
"""
function MinskerNdaoud(x, k, p)
    n = length(x)
    @assert n % k == 0 "Check that n=$n is a mutliple of k=$k"
    m = n ÷ k
    μ̄ = [mean(@view(x[(1+(l-1)*m):(l*m)])) for l in 1:k]
    σ̂p = [stdm(@view(x[(1+(l-1)*m):(l*m)]), μ̄[l], corrected = false) for l in 1:k] .^ p
    return mean(μ̄ ./ σ̂p) * harmmean(σ̂p)
end
k_mn(δ) = ceil(Int, -log(δ / 3))
struct δMinskerNdaoud{T,F} <: δ_MeanEstimator
    δ::T
    p::F
end

"""
    mean(A::AbstractArray, MeanEstimator::δMinskerNdaoud; kwargs...)

Reference: Robust and efficient mean estimation: an approach based on the properties of self-normalized sums
"""
function mean(A::AbstractArray, MeanEstimator::δMinskerNdaoud; kwargs...)
    k = k_mn(MeanEstimator.δ)
    return MinskerNdaoud(A, k, MeanEstimator.p; kwargs...)
end