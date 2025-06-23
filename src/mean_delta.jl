"""
    mean(A::AbstractArray, δ::Real, Estimator::EmpiricalMean)
    
The usual empirical mean estimator, for convenience, we authorize the extra dummy argument δ.
"""
function mean(A::AbstractArray, δ::Real, Estimator::EmpiricalMean)
    return mean(A)
end

k_mom(δ) = ceil(Int, 8log(1 / δ))
"""
    mean(A::AbstractArray, Estimator::MoM)

"""
function mean(A::AbstractArray, δ::Real, Estimator::MedianOfMean; check_multiple = false)
    k = k_mom(δ)
    if check_multiple
        n = length(A)
        if n % k != 0
            @warn "Check that n=$(n) is a mutliple of k=$k"
        end
    end
    return mean(A, MedianOfMean(k))
end

"""
    TrimmedMean(z, ε) 

"""
function TrimMean(z, ε)
    n2 = length(z)
    @assert ε < 1 "ε > 1: you cannot remove more samples than there is"
    @assert iseven(n2)
    n = n2 ÷ 2
    y = z[1:n]
    x = z[n+1:n2]
    sort!(y)
    i_min = ceil(Int, ε * n)
    i_max = n - floor(Int, ε * n)
    α = y[i_min]
    β = y[i_max]
    return TrimMean(x, α, β)
end
ε_trimmed(δ, n) = 32log(8 / δ) / (3n)

"""
     mean(A::AbstractArray, Estimator::TrimmedMean)

"""
function mean(A::AbstractArray, δ::Real, Estimator::TrimmedMean)
    ε = ε_trimmed(δ, length(A))
    return TrimMean(A, ε)
end

struct Catoni{F} <: RobustMean
    σ::F
end

α_Catoni(δ, n, σ) = sqrt((2log(2 / δ)) / (n * (1 + (2log(2 / δ)) / (n - 2log(2 / δ))))) / σ

"""
   mean(A::AbstractArray, Estimator::Catoni, kwargs...)

Reference: Catoni 
"""
function mean(A::AbstractArray, δ::Real, Estimator::Catoni, kwargs...)
    n = length(A)
    α = α_Catoni(δ, n, Estimator.σ)
    z = Z_Estimator(α, ψ_Catoni)
    return mean(A, z; kwargs...)
end

function mean(A::AbstractArray, δ::Real, Estimator::Catoni{<:Nothing}, kwargs...)
    n = length(A)
    α = α_Catoni(δ, n, std(A))
    z = Z_Estimator(α, ψ_Catoni)
    return mean(A, z; kwargs...)
end

struct Huber{F} <: RobustMean
    σ::F
end

α_Huber(δ, n, σ) = sqrt(2log(4 / δ) / n) / σ #sqrt(2 * (log(2 / δ) - log(1 - exp(-n / 2) / δ)) / n) / σ

"""
   mean(A::AbstractArray, Estimator::Huber, kwargs...)

Reference: Huber 
"""
function mean(A::AbstractArray, δ::Real, Estimator::Huber, kwargs...)
    n = length(A)
    α = α_Huber(δ, n, Estimator.σ)
    z = Z_Estimator(α, ψ_Huber)
    mean(A, z; kwargs...)
end

"""
    LeeVal(x, δ; α₀ = 0.0) 

Reference: *Optimal Sub-Gaussian Mean Estimation in ℝ* by Lee and Valiant 
"""
function LeeVal(x, δ; α₀=0.0)
    κ = MoM(x, ceil(Int, log(1 / δ)))
    f(α) = sum(min(α * (x[i] - κ)^2, 1) for i in eachindex(x)) - 1 / 3 * log(1 / δ)
    α = find_zero(f, α₀)
    return κ + mean((x[i] - κ) * (1 - min(α * (x[i] - κ)^2, 1)) for i in eachindex(x))
end
struct LeeValiant <: RobustMean end

"""
    mean(A::AbstractArray, Estimator::δLeeValiant; kwargs...)

Reference: *Optimal Sub-Gaussian Mean Estimation in ℝ* by Lee and Valiant 
"""
function mean(A::AbstractArray, δ::Real, Estimator::LeeValiant; check_multiple = false, kwargs...)
    k = ceil(Int, log(1 / δ))
    if check_multiple
        n = length(A)
        if n % k != 0
            @warn "Check that n=$(n) is a mutliple of k=$k"
        end
    end
    return LeeVal(A, δ; kwargs...)
end

"""
    MinNda(x, k::Integer, p)

Compute the Median of Minsker-Ndaoud robust mean with `k` groups of equal size (it does not permute the samples)
Leftover samples are ignored.
Reference: *Robust and efficient mean estimation: an approach based on the properties of self-normalized sums* by Minsker Ndaoud
"""
function MinNda(x, k::Integer, p)
    n = length(x)
    m = n ÷ k
    μ̄ = [mean(@view(x[(1+(l-1)*m):(l*m)])) for l in 1:k]
    σ̂p = [stdm(@view(x[(1+(l-1)*m):(l*m)]), μ̄[l], corrected=false) for l in 1:k] .^ p
    return mean(μ̄ ./ σ̂p) * harmmean(σ̂p)
end

k_mn(δ) = ceil(Int, -log(δ / 3))
struct MinskerNdaoud{T} <: RobustMean
    p::T
end

"""
    mean(A::AbstractArray, δ::Real, Estimator::MinskerNdaoud; check_multiple = false, kwargs...)

Compute the Median of Minsker Ndaoud robust mean with `k = ceil(Int, -log(δ / 3))` groups of equal size (it does not permute the samples).
Leftover samples are ignored.

`check_multiple`: Return an warning if the length of `A` is not a multiple of `k`

Reference: Robust and efficient mean estimation: an approach based on the properties of self-normalized sums
"""
function mean(A::AbstractArray, δ::Real, Estimator::MinskerNdaoud; check_multiple = false, kwargs...)
    k = k_mn(δ)
    if check_multiple
        n = length(A)
        if n % k != 0
            @warn "Check that n=$(n) is a mutliple of k=$k"
        end
    end
    return MinNda(A, k, Estimator.p; kwargs...)
end