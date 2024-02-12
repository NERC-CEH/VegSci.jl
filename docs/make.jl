using EcoVeg
using Documenter

ex_meta = quote
    # Import module(s):
    using EcoVeg

    # Data:

    # Model:
end

DocMeta.setdocmeta!(ConformalPrediction, :DocTestSetup, ex_meta; recursive=true)

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
        "Reference" => "reference.md"
    ],
)

deploydocs(;
    repo="github.com/ZekeMarshall/EcoVeg.jl.git",
    devbranch="main",
)
