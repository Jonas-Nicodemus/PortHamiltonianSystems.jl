using PortHamiltonianSystems
using Documenter, DocumenterCitations

DocMeta.setdocmeta!(PortHamiltonianSystems, :DocTestSetup, :(using PortHamiltonianSystems); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__, "..", "CITATION.bib"))

makedocs(;
    modules=[PortHamiltonianSystems],
    authors="Jonas Nicodemus <jonas.nicodemus@icloud.com> and contributors",
    repo="https://github.com/Jonas-Nicodemus/PortHamiltonianSystems.jl/blob/{commit}{path}#{line}",
    sitename="PortHamiltonianSystems.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Jonas-Nicodemus.github.io/PortHamiltonianSystems.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API" => "API.md",
    ],
    plugins=[bib],
)

deploydocs(;
    repo="github.com/Jonas-Nicodemus/PortHamiltonianSystems.jl",
    devbranch="main",
)
