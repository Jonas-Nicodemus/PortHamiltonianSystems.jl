"""
    phminreal(Σph::PortHamiltonianStateSpace; trunc_tol=1e-12)

Computes a structure preserving minimal realization of a port-Hamiltonian system [BMS22](@cite).
"""
function phminreal(Σph::PortHamiltonianStateSpace; trunc_tol=1e-12)
    # controllability Gramian X is enforced to be Q^{-1}
    Y = gram(Σph, :o)
    X = cholesky(Σph.Q) \ I

    Lx = lrcholesky(X; trunc_tol=trunc_tol)
    Ly = lrcholesky(Y; trunc_tol=trunc_tol)

    _, σ, V = svd(Lx' * Ly)
    r = length(σ[findall(σ./σ[1] .> trunc_tol)])
    V = V[:, 1:r]
    σ = σ[1:r]
    
    Wr = Ly * V * Diagonal(σ .^ (-1//2))

    return phss(Wr' * Σph.J * Wr, Wr' * Σph.R * Wr, Diagonal(inv.(σ)), Wr' * Σph.G, Wr' * Σph.P, Σph.S, Σph.N)
end