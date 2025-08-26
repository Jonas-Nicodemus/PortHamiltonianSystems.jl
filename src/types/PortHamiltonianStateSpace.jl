import Base: -
import ControlSystemsBase: StateSpace, to_matrix, AbstractNumOrArray

"""
PortHamiltonianStateSpace{T}

An object representing a port-Hamiltonian state-space system.


    dx(t)/dt = (J-R)Qx(t) + (G-P)u(t)
    y(t)     = (G+P)'Qx(t) + (S-N)u(t)


See the function [`phss`](@ref) for a user facing constructor.

# Fields:
- `J::SkewHermitian{T}`
- `R::Symmetric{T}`
- `Q::Symmetric{T}`
- `G::Matrix{T}`
- `P::Matrix{T}`
- `S::Symmetric{T}`
- `N::SkewHermitian{T}`
"""
struct PortHamiltonianStateSpace{T}
    J::SkewHermitian{T}
    R::Symmetric{T}
    Q::Symmetric{T}
    G::Matrix{T}
    P::Matrix{T}
    S::Symmetric{T}
    N::SkewHermitian{T}
end

function PortHamiltonianStateSpace(J::AbstractNumOrArray, R::AbstractNumOrArray, Q::AbstractNumOrArray, 
    G::AbstractNumOrArray, P::AbstractNumOrArray, S::AbstractNumOrArray, N::AbstractNumOrArray)
    T = promote_type(eltype(J), eltype(R), eltype(Q), eltype(G), eltype(P), eltype(S), eltype(N))

    return PortHamiltonianStateSpace{T}(
        SkewHermitian(to_matrix(T, J)), 
        Symmetric(to_matrix(T, R)), 
        Symmetric(to_matrix(T, Q)), 
        to_matrix(T, G), 
        to_matrix(T, P), 
        Symmetric(to_matrix(T, S)), 
        SkewHermitian(to_matrix(T, N)))
end

"""
    Σph = phss(J, R, Q, G, P, S, N)
    Σph = phss(J, R, Q, G)
    Σph = phss(Γ, W, Q)  

Creates a port-Hamiltonian state-space model `Σph::PortHamiltonianStateSpace{T}`
with matrix element type `T`. 
Note that `J`, `N` must be skew-symmetric and `R`, `Q`, `S` must be symmetric.
When using the structure and dissipation matrices `Γ` and `W`, they must also satisfy these properties,
i.e., `Γ` must be skew-symmetric and `W` must be symmetric.

"""
phss(args...;kwargs...) = PortHamiltonianStateSpace(args...;kwargs...)

function phss(J::AbstractNumOrArray, R::AbstractNumOrArray, Q::AbstractNumOrArray, G::AbstractNumOrArray)
    T = promote_type(eltype(J), eltype(R), eltype(Q), eltype(G))
    return phss(J, R, Q, G, zeros(T, size(G,1), size(G,2)), zeros(T, size(G,2), size(G,2)), zeros(T, size(G,2), size(G,2)))
end

function phss(Γ::AbstractMatrix, W::AbstractMatrix, m::Int)
    n = size(Γ,1) - m
    
    return phss(Γ, W, I(n))
end

function phss(Γ::AbstractMatrix, W::AbstractMatrix, Q::AbstractMatrix)
    n = size(Q,1)
    J = Γ[1:n, 1:n]
    G = Γ[1:n, n+1:end]
    N = Γ[n+1:end, n+1:end]

    R = W[1:n, 1:n]
    P = W[1:n, n+1:end]
    S = W[n+1:end, n+1:end]
    return phss(J, R, Q, G, P, S, N)
end

# Functions for number of inputs, outputs and states
ninputs(sys::PortHamiltonianStateSpace) = size(sys.G, 2)
noutputs(sys::PortHamiltonianStateSpace) = size(sys.G, 2)
nstates(sys::PortHamiltonianStateSpace) = size(sys.G, 1)

"""
    Γ(Σ) 

Returns the structure matrix of the pH system.
"""
function Γ(Σ::PortHamiltonianStateSpace)
    return SkewHermitian([Σ.J Σ.G; -Σ.G' Σ.N])
end

"""
    W(Σ)

Returns the dissipation matrix of the pH system.
"""
function W(Σ::PortHamiltonianStateSpace)
    return Symmetric([Σ.R Σ.P; Σ.P' Σ.S])
end

function Base.getproperty(Σ::PortHamiltonianStateSpace, s::Symbol)
    if s === :nx
        return nstates(Σ)
    elseif s === :nu
        return ninputs(Σ)
    elseif s === :ny
        return noutputs(Σ)
    elseif s === :Γ
        return Γ(Σ)
    elseif s === :W
        return W(Σ)
    else
        return getfield(Σ, s)
    end
end

function -(sys1::PortHamiltonianStateSpace, sys2::PortHamiltonianStateSpace)
    return ss(sys1) - ss(sys2)
end

function -(sys1::StateSpace, sys2::PortHamiltonianStateSpace)
    return sys1 - ss(sys2)
end

function -(sys1::PortHamiltonianStateSpace, sys2::StateSpace)
    return ss(sys1) - sys2
end

_string_mat_with_headers(X) = _string_mat_with_headers(Matrix(X))
function _string_mat_with_headers(X::Matrix)
    p = (io, m) -> Base.print_matrix(io, m)
    return replace(sprint(p, X), "\"" => " ")
end

Base.print(io::IO, sys::PortHamiltonianStateSpace) = show(io, sys)

function Base.show(io::IO, sys::PortHamiltonianStateSpace)
    println(io, typeof(sys))
    println(io, "J = \n", _string_mat_with_headers(sys.J))
    println(io, "R = \n", _string_mat_with_headers(sys.R))
    println(io, "Q = \n", _string_mat_with_headers(sys.Q))
    println(io, "G = \n", _string_mat_with_headers(sys.G))
    println(io, "P = \n", _string_mat_with_headers(sys.P))
    println(io, "S = \n", _string_mat_with_headers(sys.S))
    println(io, "N = \n", _string_mat_with_headers(sys.N))
end