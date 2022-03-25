using RobustMeans
using Test

@testset "RobustMean.jl" begin
using Distributions
using Statistics, Random
Random.seed!(3)
n = 3 * 7 * 8
M = 10^5
b = 2.0001
dist = Pareto(b)
μ = mean(dist)
σ = std(dist)
x = rand(dist, M, n)

δ = 3exp(-8)
p = 1 # for Minsker Ndaoud
estimator = [EmpiricalMean(), MedianOfMean(1), TrimmedMean(0, 0), Catoni(σ), Huber(σ), LeeValiant(), MinskerNdaoud(p)]

μ̂ = zeros(M, length(estimator))
for (i, es) in enumerate(estimator)
    for m in 1:M
        μ̂[m, i] = mean(x[m, :], δ, es)
    end
end

for (i, es) in enumerate(estimator)
    if typeof(es) <: MinskerNdaoud
        continue # constant in the bound not explciit.
    else
        @test quantile(μ̂[:, i], 1 - δ) < σ * bound(δ, n, es)
    end
end

end
