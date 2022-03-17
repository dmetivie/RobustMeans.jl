using RobustMean
using Test

@testset "RobustMean.jl" begin
    using Distribution

    n = 40
    M = 10^6
    b = 3.0001
    dist = Pareto(b)
    μ = mean(dist)
    σ = std(dist)
    x = rand(dist, M, n)

    
    @test isapprox(β, probs(mix_mle)[1]; rtol = rtol)
    @test isapprox(θ₁, p[1]...; rtol = rtol)
    @test isapprox(α, p[2][1]; rtol = rtol)
    @test isapprox(θ₂, p[2][2]; rtol = rtol)
end
