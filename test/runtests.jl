using PortHamiltonianSystems
using Test, Aqua

@testset "PortHamiltonianSystems.jl" begin
    Aqua.test_all(PortHamiltonianSystems; piracies=false)

    include("test_types.jl")
    include("test_analysis.jl")
    include("test_convert.jl")
    include("test_gramians.jl")
    include("test_freqresp.jl")
    include("test_kyp.jl")
    include("test_linalg.jl")
    include("test_minimality.jl")
end
