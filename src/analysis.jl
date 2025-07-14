import LinearAlgebra: norm

"""
    norm(Σph, p=2; kwargs...)

Converts a `PortHamiltonianStateSpace` to a `StateSpace` and calls `norm` on it.
For more details see [ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl) package.
"""
function norm(sys::PortHamiltonianStateSpace, p::Real=2; kwargs...)
    return norm(ss(sys), p; kwargs...)
end

"""
    ispassive(Σ::StateSpace; opt=:lmi, kwargs...)
    ispassive(Σ::PortHamiltonianStateSpace; kwargs...)
    
Checks whether the system `Σ` is passive by solving the KYP inequality using [`kyp`](@ref) (if opt=`:lmi`) 
or checking the Popov function for passivity violations via [`sampopov`](@ref) (if opt=`:popov`). 
"""
function ispassive(Σ::StateSpace; opt=:lmi, kwargs...)
    if opt==:lmi
        return ispassive_lmi(Σ; kwargs...)
    elseif opt==:popov
        return ispassvie_sampopov(Σ; kwargs...)
    else 
        error("opt must be either :lmi or :popov")
    end
end


function ispassive(Σ::PortHamiltonianStateSpace; kwargs...)
    return ispassive(ss(Σ); kwargs...)
    
end

"""

    sampopov(Σ; ω=10 .^ range(-15, stop=5, length=5000))
    
Samples the Popov function for the system `Σ` at ranges of frequencies `ω`.
"""
function sampopov(Σ; ω=10 .^ range(-15, stop=5, length=5000))
    Φ = popov(Σ, ω)
    F = eigen.([Φ[:, :, i] for i ∈ axes(Φ, 3)])
    Λmin = [minimum(real.(F[i].values)) for i ∈ eachindex(F)]
    
    c = findall(i -> (sign(Λmin[i]) != sign(Λmin[i+1])), range(1, length(Λmin)-1))
    c = [0; c...; length(Λmin)]
    Ωi = [c[i]+1:c[i+1] for i ∈ range(1, length(c)-1)]
    
    if Λmin[1] > 0
        Ωp = Ωi[1:2:end]
        Ωnp = Ωi[2:2:end] 
    else
        Ωp = Ωi[2:2:end]
        Ωnp = Ωi[1:2:end]
    end

    return Λmin, ω, Ωp, Ωnp, F, c
end

function ispassvie_sampopov(Σ; kwargs...)
    _, _, _, Ωnp, _, _  = sampopov(Σ; kwargs...)
    return isempty(Ωnp)
end

function ispassive_lmi(Σ::StateSpace; kwargs...)
    X = kyp(Σ; kwargs...)
    return ispsd(X) && ispsd(kypmat(Σ, X))
end
