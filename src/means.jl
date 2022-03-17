"""
    MoM(x::AbstractArray, k::Integer)
    Compute the Median of Mean with `k` groups (it does not permute the samples)
"""
function MedianOfMean(x::AbstractArray, k::Integer)
    n = length(x)
    @assert n % k == 0 "Check that n=$n is a mutliple of k=$k"
    m = n ÷ k
    return median(mean(@view(x[(1+(l-1)*m):(l*m)])) for l = 1:k)
end

"""
    TrimmedMean(x::AbstractArray, k::Integer)
    Compute the Trimmed Mean thresolding data smaller `α` or larger than `β`
"""
function φ(x, α, β)
    if x < α
        return α
    elseif x > β
        return β
    else
        x
    end
end

TrimmedMean(x::AbstractArray, α, β) = mean(φ(x[i], α, β) for i in eachindex(x))


"""
 Z-Estimators
"""
ψ_Catoni(x) = x < 0 ? -log(1 - x + x^2 / 2) : log(1 + x + x^2 / 2)
ψ_Huber(x) = abs(x) ≤ 1 ? x : sign(x)

Rα(y, x::AbstractArray, α, ψ::Function) = sum(ψ(α * (x[i] - y)) for i in eachindex(x))

function Z_estimator(x::AbstractArray, α, ψ::Function; ini = median(x))
    #! For some reason it seems that `find_zero` only accepts `Float64`
    f(y) = Rα(y, x, α, ψ)
    return find_zero(f, ini)
end

