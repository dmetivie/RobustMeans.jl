using RobustMean
using Test

@testset "RobustMean.jl" begin
    using Distributions
    using Statistics, Random
    Random.seed!(3)
    n = 3 * 7 * 8
    M = 10^6
    b = 2.0001
    dist = Pareto(b)
    μ = mean(dist)
    σ = std(dist)
    x = rand(dist, M, n)

    δ = 3exp(-8)
    p = 1 # for Minsker Ndaoud
    estimator = [EmpiricalMean(), δMoM(δ), δTrimM(δ), δCatoni(δ, σ), δHuber(δ, σ), δLeeValiant(δ), δMinskerNdaoud(δ, p)]

    μ̂ = zeros(M, length(estimator))
    for (i, es) in enumerate(estimator)
        for m in 1:M
            μ̂[m, i] = mean(x[m, :], es)
        end
    end

    for (i, es) in enumerate(estimator)
        if typeof(es) == δMinskerNdaoud
            continue # constant in the bound not explciit.
        else
            @test quantile(μ̂[:, i], 1 - δ) < σ * bound(n, δ, es)
        end
    end

end
