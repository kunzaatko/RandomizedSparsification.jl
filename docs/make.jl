using RandomizedSparsification
using Documenter, DocumenterCitations

DocMeta.setdocmeta!(
    RandomizedSparsification,
    :DocTestSetup,
    :(using RandomizedSparsification);
    recursive=true,
)

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    # style=:authoryear
)

makedocs(;
    modules=[RandomizedSparsification],
    authors="Martin Kunz <martinkunz@email.cz> and contributors",
    sitename="RandomizedSparsification.jl",
    format=Documenter.HTML(;
        canonical="https://kunzaatko.github.io/RandomizedSparsification.jl",
        edit_link="trunk",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
    plugins=[bib],
)

deploydocs(; repo="github.com/kunzaatko/RandomizedSparsification.jl", devbranch="trunk")
