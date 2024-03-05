using VegSci
using Documenter

ex_meta = quote
    # Import module(s):
    using VegSci

    # Data:

    # Model:
end

DocMeta.setdocmeta!(VegSci, :DocTestSetup, ex_meta; recursive=true)

makedocs(;
    modules=[VegSci],
    authors="Zeke Marshall",
    sitename="VegSci.jl",
    format=Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical="https://ZekeMarshall.github.io/VegSci.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Methods" => [
            "Correspondence Analysis" => "methods/CorrespondenceAnalysis.md"
            "Species Fidelity" => "methods/SpeciesFidelity.md"
            # "Binary Similarity" => "BinarySimilarity.md"
        ],
        "Reference" => "reference.md"
    ],
)

deploydocs(;
    repo="github.com/ZekeMarshall/VegSci.jl.git",
    devbranch="main",
)
