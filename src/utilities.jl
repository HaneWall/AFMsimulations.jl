
""" 
returns effective Young's module of two interacting materials
"""
function eff_young_module(s::Sample, t::Cantilever)
    E = ( (1-s.ν^2)/s.E  + (1-t.ν^2)/t.E )^(-1)
    return E
end


"""
returns force-distance relation 
"""
function force_distance(x::Float64, exp::AFM_DMT_experiment)
    R = exp.tip.R
    H = exp.sample.H
    E = exp.eff_young
    d = exp.d
    a_0 = exp.a_0
    if x < d - a_0
        return -H*R/(6*(d-x)^2)
    else 
        return -H*R/(6 * a_0^2) + 4/3*E*sqrt(R)* (x - (d - a_0))^(3/2) 
    end
end

"""
implement lock in fast fourier trasnform --> get phase and amplitude 
"""