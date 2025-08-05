# PortHamiltonianSystems.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Jonas-Nicodemus.github.io/PortHamiltonianSystems.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Jonas-Nicodemus.github.io/PortHamiltonianSystems.jl/dev/)
[![Build Status](https://github.com/Jonas-Nicodemus/PortHamiltonianSystems.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Jonas-Nicodemus/PortHamiltonianSystems.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/jonas-nicodemus/PortHamiltonianSystems.jl/graph/badge.svg?token=7CZQ640M8K)](https://codecov.io/gh/jonas-nicodemus/PortHamiltonianSystems.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A package for port-Hamiltonian systems in Julia.

The considered systems are of the form

```math
    \Sigma\quad\left\{\quad\begin{aligned}
        \dot{x}(t) &= (J - R)Qx(t) + (G - P)u(t), \\
        y(t) &= (G + P)^\top Qx(t) + (S - N)u(t)
    \end{aligned}\right.,
```
where $J=-J^\top\in \mathbb{R}^{n\times n}$ is skew-symmetric, $R\in \mathbb{R}^{n\times n}$ is positive semi-definite, $Q\in \mathbb{R}^{n\times n}$ is positive definite, $G\in \mathbb{R}^{n\times m}$, $P\in \mathbb{R}^{n\times m}$, $S\in \mathbb{R}^{m\times m}$ and $N\in \mathbb{R}^{m\times m}$.

## Installation

Install with the Julia package manager [Pkg](https://pkgdocs.julialang.org/):
```julia
pkg> add PortHamiltonianSystems # Press ']' to enter the Pkg REPL mode.
```
or
```julia
julia> using Pkg; Pkg.add("PortHamiltonianSystems")
```

## Documentation

Some available commands are:
##### Constructing systems
`phss`
##### Analysis
`norm, ispassive, sampopov`
##### Gramians
`grampd, gram, prgrampd, prgram`
##### Minimal realization
`phminreal`
##### Frequency response
`popov`

### Example

```julia
using LinearAlgebra, ControlSystemsBase
using PortHamiltonianSystems

J = [0 -1; 1 0]
R = [1 -1; -1 2]
Q = I(2)
G = [1; 0;;]

Σph = phss(J, R, Q, G)
# This generates the system:
# PortHamiltonianStateSpace{Float64}
# J = 
#   0.0  -1.0
#  1.0  0.0
# R = 
#   1.0  -1.0
#  -1.0   2.0
# Q = 
#  1.0  0.0
#  0.0  1.0
# G = 
#  1.0
#  0.0
# P = 
#  0.0
#  0.0
# S = 
#  0.0
# N = 
#  0.0

# Compute the observability Gramian
gram(Σph, :o)
# 2×2 Matrix{Float64}:
#  0.5  0.0
#  0.0  0.0

# Compute a numerically minimal realization
Σphr = phminreal(Σph)

# Compute H2 and Hinf errors
norm(Σph - Σphr) # 9.56908750203461e-17
norm(Σph - Σphr, Inf) # 1.8182097724708874e-16

# Computes unstructured state-space realization
Σ = ss(Σph)

# Checking passivity
ispassive(Σ) # true

# Computes pH realization for a given posdef X satisfying the KYP inequality 
X = kyp(Σ)
Σph2 = phss(Σ, X)

norm(Σph - Σph2) # 1.3470909371896273e-16
```
