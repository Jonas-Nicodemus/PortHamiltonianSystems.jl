module Test_analysis

using Test

using LinearAlgebra, ControlSystemsBase
using PortHamiltonianSystems

@testset "test_analysis.jl" begin
    A = [-2 1.; -1. -1.]
    B = [6.; 0]
    C = B'
    D = 1
    Σp = ss(A,B,C,D)

    A = [-0.08 0.83 0 0;-0.83 -0.08 0 0;0 0 -0.7 9;0 0 -9 -0.7]
    B = [1 1;0 0;1 -1;0 0]
    C = [0.4 0.6;0 0;0.4 1;0 0]' 
    D = [0.3 0;0 -0.15]
    Σnp = ss(A, B, C, D)

    @testset "test_ispassive" begin
        # for opt ∈ [:lmi, :popov, :ham]
        for opt ∈ [:lmi, :popov]
            @test ispassive(Σp; opt=opt) == true
            @test ispassive(Σnp; opt=opt) == false
        end
    end
end

end