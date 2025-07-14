module Test_minimality

using Test

using PortHamiltonianSystems

using LinearAlgebra, ControlSystemsBase

@testset "test_minimality.jl" begin

    A = [-1 0; 2 -2]
    B = [1; 0]
    Q = Array(I(2))
    Σ = ss(A, B, B', 0)

    Σph = phss(Σ, Q)
   
    @testset "phminreal" begin
        Σm = phminreal(Σph)
        r = rank(gram(Σ, :o))
        @test size(Σm.J) == (r, r)
        @test norm(Σph - Σm) < 1e-12
    end
end

end