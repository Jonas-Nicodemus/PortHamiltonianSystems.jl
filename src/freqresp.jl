"""
    popov(Σ, s)
    popov(Σ, ω)

Evaluates the popov function of the system `Σ` at the complex variable `s`.

    Φ(s) = Σ(s) + Σ(-s)'

where `Σ(s)` is the frequency response (transfer function) of `Σ` at the complex variable `s`.
"""
function popov(Σ::StateSpace, s::Complex)
    return evalfr(Σ, s) + transpose(evalfr(Σ, -s))
end

function popov(Σ::StateSpace, ω::AbstractVector{W}) where W <: Real
    return popov(freqresp(Σ, ω))
end

function popov(Σ::PortHamiltonianStateSpace, s)
    return popov(ss(Σ), s)
end

function popov(H)
    return real(H + permutedims(conj(H), [2,1,3]))
end