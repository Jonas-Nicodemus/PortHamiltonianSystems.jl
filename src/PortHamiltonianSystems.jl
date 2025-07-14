module PortHamiltonianSystems

using LinearAlgebra, SkewLinearAlgebra, VectorizationTransformations
using ControlSystemsBase, MatrixEquations
using JuMP, Hypatia

export PortHamiltonianStateSpace, phss # types/PortHamiltonianStateSpace.jl
export ispassive, sampopov # analysis.jl
export grampd, gram, prgrampd, prgram # gramians.jl
export kyp, kypare, kypmat, kypmin, kypmax # kyp.jl
export compose, decompose # convert.jl
export popov # freqresp.jl
export project_psd, ispsd, lrcholesky, tsvd # linalg.jl
export unvech, unvec # linalg.jl
export phminreal # minimality.jl

include("types/PortHamiltonianStateSpace.jl")
include("analysis.jl")
include("convert.jl")
include("freqresp.jl")
include("gramians.jl")
include("kyp.jl")
include("linalg.jl")
include("minimality.jl")

end
