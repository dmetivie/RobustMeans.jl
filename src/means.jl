"""
    mean(A::AbstractArray, Estimator::EmpiricalMean)
    
The usual empirical mean estimator. Nothing fancy.
"""
struct EmpiricalMean <: MeanEstimator end
function mean(A::AbstractArray, Estimator::EmpiricalMean)
    return mean(A)
end

"""
MedianOfMean(x::AbstractArray, k::Integer)

Compute the Median of Mean with `k` groups (it does not permute the samples)
"""
function MoM(x::AbstractArray, k::Integer)
    n = length(x)
    @assert n % k == 0 "Check that n=$n is a mutliple of k=$k"
    m = n ÷ k
    return median(mean(@view(x[(1+(l-1)*m):(l*m)])) for l = 1:k)
end
struct MedianOfMean{T<:Integer} <: RobustMean
    k::T # number of blocks
end
"""
mean(A::AbstractArray, Estimator::MoM)

The Median of Mean estimator.
"""
function mean(A::AbstractArray, Estimator::MedianOfMean)
    return MoM(A, Estimator.k)
end

struct TrimmedMean{T<:Integer} <: RobustMean
    α::T # 
    β::T
end
"""
    TrimmedMean(x::AbstractArray, k::Integer)

Compute the Trimmed Mean thresolding data smaller `α` or larger than `β`
"""
TrimMean(x::AbstractArray, α, β) = mean(φ(x[i], α, β) for i in eachindex(x))
function φ(x, α, β)
    if x < α
        return α
    elseif x > β
        return β
    else
        x
    end
end
"""
mean(A::AbstractArray, Estimator::TrimmedMean)

The Trimmed Mean estimator. (p)
"""
function mean(A::AbstractArray, Estimator::TrimmedMean)
    return TrimMean(A, Estimator.α, Estimator.β)
end

"""
    Z_estimator(x::AbstractArray, α, ψ::Function; ini = median(x))

Implement Z estimators given a 
"""
function Z_estimator(x::AbstractArray, α, ψ::Function; ini=median(x))
    #! For some reason it seems that `find_zero` only accepts `Float64`
    f(y) = Rα(y, x, α, ψ)
    return find_zero(f, ini)
end

Rα(y, x::AbstractArray, α, ψ::Function) = sum(ψ(α * (x[i] - y)) for i in eachindex(x))

struct Z_Estimator{F,T} <: MeanEstimator
    α::T # scaling parameter
    ψ::F # Influence function
end
"""
mean(A::AbstractArray, Estimator::Z_Estimator)

A Z estimator given an influence function `x->ψ(x)` and a scaling parameter `α`.
"""
function mean(A::AbstractArray, Estimator::Z_Estimator; kwargs...)
    return Z_estimator(A, Estimator.α, Estimator.ψ; kwargs...)
end

"""
    ψ_Catoni(x)

Catoni's influence function
"""
ψ_Catoni(x) = x < 0 ? -log(1 - x + x^2 / 2) : log(1 + x + x^2 / 2)

"""
    ψ_Huber(x)

ψ_Huber's influence function
"""
ψ_Huber(x) = abs(x) ≤ 1 ? x : sign(x)