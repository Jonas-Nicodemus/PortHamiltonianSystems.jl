module Test_freqresp

using Test

using PortHamiltonianSystems

using LinearAlgebra, ControlSystemsBase

@testset "test_freqresp.jl" begin

    A = [-2. 1.; -1. -1.]
    B = [6.; 0;]
    C = B'
    D = [1.;;]

    Σ = ss(A,B,C,D)
        
    @testset "popov" begin
        @testset "s" begin
            s = 1 + 1im
            @test norm(popov(Σ, s) - (evalfr(Σ, s) + transpose(evalfr(Σ, -s)))) < 1e-8
        end

        @testset "ω" begin
            ω = 10 .^ range(-2, stop=2, length=5)
            @test norm(popov(Σ, ω) - (freqresp(Σ, ω) + permutedims(freqresp(Σ, -ω), [2,1,3]))) < 1e-8
        end            
    end
end

end