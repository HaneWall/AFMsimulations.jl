
""" 
returns effective Young's module of two interacting materials
"""
function eff_young_module(s::Sample, t::Cantilever)
    E = ( (1-s.ν^2)/s.E  + (1-t.ν^2)/t.E )^(-1)
    return E
end