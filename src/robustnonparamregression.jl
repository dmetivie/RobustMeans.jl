function r̂ₙ(x, X::AbstractArray, Y::AbstractArray, kernel::SmoothingKernel)
    Nₕ = [kernel.f((X[i] - x) / kernel.h) for i in eachindex(X)]
    sumNₕ = sum(Nₕ)
    if sumNₕ == 0
        return 0
    else
        return sum(Nₕ[i] * Y[i] for i in eachindex(Y)) / sum(Nₕ)
    end
end

function RobustMeans.r̂ₙ(x::AbstractArray, X::AbstractArray, Y::AbstractArray, kernel::SmoothingKernel)
	out = zeros(eltype(Y), length(x))
	Nₕ = similar(Y)
	for (j, xx) in enumerate(x)
    	Nₕ[:] = [kernel.f((X[i] - xx) / kernel.h) for i in eachindex(X)]
    	sumNₕ = sum(Nₕ)
    	if sumNₕ ≈ 0
        	continue
    	else
        	out[j] = sum(Nₕ[i] * Y[i] for i in eachindex(Y)) / sum(Nₕ)
    	end
	end
	return out
end

function r̂ₙMOM(x, X, Y, kernel, m)
    N = length(X)
    mom = zeros(m)
    for (i, c) in enumerate(RobustMeans.chunk(N, m))
        mom[i] = r̂ₙ(x, X[c], Y[c], kernel)
    end
    return median(mom)
end