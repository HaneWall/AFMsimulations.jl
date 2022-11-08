""" 
R .. radius of tip, Q .. quality factor of tip, k .. spring constant, 
E .. Young's module, ν .. poisson coefficient, f_0 .. eigenfrequency. 
"""
struct Cantilever
    R::Float64
    Q::Float64
    k::Float64
    E::Float64
    ν::Float64
    f_0::Float64
end


""" 
H .. Hamaker constant, E .. Young's module, ν .. poisson coefficient.
"""
struct Sample
    H::Float64
    E::Float64
    ν::Float64
end

"""
Lennard Jones experiment
"""
struct AFM_LJ_experiment
    sample::Sample
    tip::Cantilever
    σ::Float64
    V_0::Float64
    d::Float64
    function AFM_LJ_experiment(sam::Sample, tip::Cantilever, σ::Float64, V_0::Float64, d::Float64)
        new(sam, tip, σ, V_0, d) 
    end
end



"""
DMT experiment, respects elastic behavior within the hertzian contact model 
"""
struct AFM_DMT_experiment
    sample::Sample
    tip::Cantilever
    a_0::Float64
    d::Float64
    eff_young::Float64
    function AFM_DMT_experiment(sam::Sample, tip::Cantilever, a_0::Float64, d::Float64)
        eff_young = eff_young_module(sam, tip)
        new(sam, tip,  a_0, d, eff_young) 
    end
end


"""
DMT experiment with exponential damping term --> eDMT
"""
struct AFM_eDMT_experiment
    sample::Sample
    tip::Cantilever
    γ_0::Float64
    x_γ::Float64
    a_0::Float64
    d::Float64
    eff_young::Float64
    function AFM_eDMT_experiment(sam::Sample, tip::Cantilever, γ_0::Float64, x_γ::Float64, a_0::Float64, d::Float64)
        eff_young = eff_young_module(sam, tip)
        new(sam, tip, γ_0, x_γ, a_0, d, eff_young) 
    end
end

"""
DMT experiment with viscious damping term --> vDMT
"""
struct AFM_vDMT_experiment
    sample::Sample
    tip::Cantilever
    η::Float64
    a_0::Float64
    d::Float64
    eff_young::Float64
    function AFM_vDMT_experiment(sam::Sample, tip::Cantilever, η::Float64, a_0::Float64, d::Float64)
        eff_young = eff_young_module(sam, tip)
        new(sam, tip, η, a_0, d, eff_young) 
    end
end