using Aqua

@testset "Aqua.jl" begin
    Aqua.test_all(
      EcoVeg;
      piracies=false,
    )
  end