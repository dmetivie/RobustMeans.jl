abstract type δ_MeanEstimator <: MeanEstimator end

k_mom(δ) = ceil(Int, 8log(1 / δ))
struct δMoM{T} <: δ_MeanEstimator
    δ::T
end
function mean(A::AbstractArray, MeanEstimator::δMoM)
    k = k_mom(δMoM.δ)
    return mean(A, MoM(k))
end

ϵ_trimmed(δ, n) = 32log(8 / δ) / (3n)
