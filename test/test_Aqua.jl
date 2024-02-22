using Aqua

@testset "Aqua.jl" begin
    Aqua.test_all(
      VegSci;
      piracies=false,
    )
  end