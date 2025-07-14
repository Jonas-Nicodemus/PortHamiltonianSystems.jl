module Test_utils

using Test

using PortHamiltonianSystems

using LinearAlgebra

@testset "test_linalg.jl" begin

    @testset "lrcholesky" begin
        B = [1 2 3 4; 5 6 7 8]; A = B'B;
        L = lrcholesky(A)
        @test size(L) == (size(A, 1), rank(A))
        @test norm(L*L' - A) < 1e-12
    end

    @testset "tsvd" begin
        B = [1 2 3 4; 5 6 7 8]; A = B'B;
        F1 = tsvd(A)
        @test size(F1.U) == (size(A, 1), rank(A))
        @test size(F1.Vt) == (rank(A), size(A, 2))
        @test length(F1.S) == rank(A)
        @test norm(Matrix(F1) - A) < 1e-12 

        F2 = tsvd(A, 2)
        @test size(F2.U) == (size(A, 1), rank(A))
        @test size(F2.Vt) == (rank(A), size(A, 2))
        @test length(F2.S) == rank(A)
        @test norm(Matrix(F2) - A) < 1e-12 
    end

    @testset "project_psd" begin
        eigtol = 1e-12
        M = [2 1; 1 -1]
        Mpsd = project_psd(M; eigtol=eigtol)
        @test norm(Mpsd - Mpsd') == 0
        @test all(eigvals(Mpsd) .>= -eigtol)
    end

    @testset "unvec" begin
        M = [1 2; 3 4] 
        v1 = vec(M)
        @test unvec(v1) == M 
        v2 = reshape(v1, 1, 4)
        @test unvec(v2) == M
    end

    @testset "unvech" begin
        s = [1, 2, 3]
        S = unvech(s)
        @test S == [1 2; 2 3]

        s = [1, 2, 3, 4, 5, 6]
        S = unvech(s)
        @test S == [1 2 3; 2 4 5; 3 5 6]
    end

end

end