using RandomizedSparsification
using Documenter

DocMeta.setdocmeta!(RandomizedSparsification, :DocTestSetup, :(using RandomizedSparsification); recursive=true)

makedocs(;
    modules=[RandomizedSparsification],
    authors="Martin Kunz <martinkunz@email.cz> and contributors",
    sitename="RandomizedSparsification.jl",
    format=Documenter.HTML(;
        canonical="https://kunzaatko.github.io/RandomizedSparsification.jl",
        edit_link="trunk",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kunzaatko/RandomizedSparsification.jl",
    devbranch="trunk",
)
