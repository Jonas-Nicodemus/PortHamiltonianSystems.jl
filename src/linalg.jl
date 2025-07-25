"""
    L = lrcholesky(X; trunc_tol=1e-12)

Computes a low-rank approximate Cholesky-like factorization of a symmetric positive semi-definite matrix ``X``
s.t. ``X = L * L'`` (up to a prescribed tolerance `trunc_tol`).
"""
function lrcholesky(X::Matrix; kwargs...)
    return lrcholesky(hermitianpart(X); kwargs...)
end
function lrcholesky(X::Hermitian; trunc_tol=1e-12)
    d,L = eigen(X)
    # remove negative eigenvalues (numerical errors)
    idx = findall(v -> v >= 0, d)
    dr, Lr = truncation(d[idx], L[:, idx]; trunc_tol)

    return Lr*diagm(sqrt.(dr))
end

"""
    truncation(d, L, trunc_tol) -> (dr, Lr)

Computes a rank revealing factorization for a given LDL-decomposition of
``S = L * \\mathrm{diag}(d) * L^T`` (up to a prescribed tolerance `trunc_tol`)
such that ``L_r * diag(d_r) * L_r^T \\approx S``.
"""
function truncation(d, L; trunc_tol=1e-12)
    Q,R = qr(L)
    tmp = Symmetric(R*diagm(d)*R')
    d,U = eigen(tmp)
    p = sortperm(d, by=abs, rev=true)
    d = d[p]
    trunc_index = findlast(abs.(d/d[1]) .>= trunc_tol)
    return d[1:trunc_index], Q*U[:,p[1:trunc_index]]
end

"""
    Mpsd = project_psd(M; eigtol=1e-8)

Returns the nearest positive semi-definite matrix to `M` by setting negative eigenvalues to `eigtol`.
"""    
function project_psd(M; eigtol=0)
    Λ, U = eigen(hermitianpart(M))
    return hermitianpart(U * diagm(clamp.(Λ, eigtol, Inf)) * U')
end

"""
    ispsd(M)

Returns `true` if `M` is positive semi-definite, otherwise `false`.
"""
function ispsd(M; ϵ=1e-8)
    return isposdef(M + ϵ*I) 
end

"""
    M = unvech(v)

Returns the Hermitian matrix `M` from the half-vectorized `v`, i.e., the inverse of `VectorizationTransforms.vech`.
"""
function unvech(v::Vector)
    n = Int((sqrt(8 * length(v) + 1) - 1) / 2)
    D = duplication_matrix(n)
    return hermitianpart(reshape(D * v, n, n))
end

"""
    M = unvec(v)

Returns the matrix `M` from the vectorized `v`, i.e., the inverse of `LinearAlgebra.vec`.
"""
function unvec(v::Vector)
    n = Int(sqrt(length(v)))
    return reshape(v, n, n)
end

function unvec(v::AbstractMatrix)
    @assert size(v, 1) == 1
    return unvec(vec(v))
end

"""
    F = tsvd(A; kwargs...)

Returns the truncated singular value decomposition of `A` by truncating small singular values below `ε`.
"""
function tsvd(A::AbstractMatrix; ε=1e-12, kwargs...)
    return tsvd(svd(A; kwargs...); ε=ε)
end
function tsvd(A::AbstractMatrix, r::Int; kwargs...)
    return tsvd(svd(A; kwargs...), r)
end

function tsvd(F::SVD; ε=1e-12)
    mask = F.S./F.S[1] .> ε
    return SVD(F.U[:, mask], F.S[mask], F.Vt[mask, :])
end

function tsvd(F::SVD, r::Int)
    return SVD(F.U[:, 1:r], F.S[1:r], F.Vt[1:r, :])
end
