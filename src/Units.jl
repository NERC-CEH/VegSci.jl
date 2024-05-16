using NamedArrays

function prop_to_domin(na::NamedArray)

  function _prop_to_domin(x::Float64)

    if x .== 0.0
      domin = 0
    elseif x .< 0.01
      domin = 1
    elseif x .< 0.025
      domin = 2
    elseif x .< 0.04
      domin = 3
    elseif x .< 0.1
      domin = 4
    elseif x .< 0.25
      domin = 5
    elseif x .< 0.33
      domin = 6
    elseif x .< 0.5
      domin = 7
    elseif x .< 0.75
      domin = 8
    elseif x .< 0.9
      domin = 9
    elseif x .< 1
      domin = 10
    end
  
    return domin

  end

  replace!(x -> _prop_to_domin(x), na)

  return na

end

function perc_to_domin(na::NamedArray)

  function _perc_to_domin(x::Float64)
    
    if x .== 0.0
      domin = 0
    elseif x .< 1.0
      domin = 1
    elseif x .< 2.5
      domin = 2
    elseif x .< 4.0
      domin = 3
    elseif x .< 10.0
      domin = 4
    elseif x .< 25.0
      domin = 5
    elseif x .< 33.0
      domin = 6
    elseif x .< 50.0
      domin = 7
    elseif x .< 75.0
      domin = 8
    elseif x .< 90.0
      domin = 9
    elseif x .< 100.0
      domin = 10
    end
  
    return domin

  end

  replace!(x -> _perc_to_domin(x), na)

  return na

end

function relfreq_to_constancy(na::NamedArray)

  function _relfreq_to_constancy(x::Float64)

    if x .< 20
      constancy = 1
    elseif x .< 40
      constancy = 2
    elseif x .< 60
      constancy = 3
    elseif x .< 80
      constancy = 4
    elseif x .< 100
      constancy = 5
    end

    return constancy

  end

  replace!(x -> _relfreq_to_constancy(x), na)

  return na

end