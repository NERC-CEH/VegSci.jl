setup_docs = quote

    # Environment:
    using Pkg
    Pkg.activate("docs")

    # Dependencies:
    using DataFrames
    using Documenter
    using LinearAlgebra
    using NamedArrays
    using PrettyTables
    using Statistics
    using Suppressor

    # Setup:
    theme(:wong)
    Random.seed!(2023)
    www_path = "$(pwd())/docs/src/www"
end;