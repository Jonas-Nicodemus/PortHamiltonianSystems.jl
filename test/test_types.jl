using Test

using LinearAlgebra
using PortHamiltonianSystems

@testset "PortHamiltonianStateSpace" begin
    J = [0. 1.; -1. 0.]
    R = [2. 0.; 0. 1.]
    Q = [1. 0.; 0. 1.]
    G = [6.; 0.;;]
    P = zero(G)
    S = [1.;;]
    N = zero(S)

    Γ = [J G; -G' N]
    W = [R P; P' S]

    @testset "phss" begin
        Σ1 = phss(J, R, Q, G, P, S, N)
        Σ2 = phss(J, R, Q, G)
        Σ3 = phss(Γ, W, Q)
        @test norm(Σ1 - Σ3) < 1e-12


        S = 1
        Σ4 = phss(J, R, Q, G, P, S, N)

        Σ5 = phss(Γ, W, 1)
        @test norm(Σ3 - Σ5) < 1e-12
    end

    @testset "Γ & W" begin
        Σ = phss(J, R, Q, G, P, S, N)
        @test norm(Γ - PortHamiltonianSystems.Γ(Σ)) < 1e-12
        @test norm(W - PortHamiltonianSystems.W(Σ)) < 1e-12
    end

    @testset "getproperty" begin
        Σ = phss(J, R, Q, G, P, S, N)

        @test Σ.nx == size(G, 1)
        @test Σ.nu == size(G, 2)
        @test Σ.ny == size(G, 2)

        @test Σ.Γ == Γ
        @test Σ.W == W
    end
end
