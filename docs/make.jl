using EcoVeg
using Documenter

DocMeta.setdocmeta!(EcoVeg, :DocTestSetup, :(using EcoVeg); recursive=true)

makedocs(;
    modules=[EcoVeg],
    authors="Zeke Marshall",
    sitename="EcoVeg.jl",
    format=Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical="https://ZekeMarshall.github.io/EcoVeg.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ZekeMarshall/EcoVeg.jl.git",
    devbranch="main",
)
