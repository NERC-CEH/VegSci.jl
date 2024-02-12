setup_docs = quote

    # Environment:
    using Pkg
    Pkg.activate("docs")

    # Dependencies:
    using EcoVeg
    using DataFrames
    using Documenter
    using LinearAlgebra
    using NamedArrays
    using PrettyTables
    using Statistics
    using Suppressor
    using Skipper

    # Setup:
    Random.seed!(2023)
    www_path = "$(pwd())/docs/src/www"
end;