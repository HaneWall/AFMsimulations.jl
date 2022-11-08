"""
Lennard-Jones model, smooth interaction potential
"""
function f_LJ!(dx, x, p, t)
    σ, Q, Ω, Γ, d, k_c, ϕ = p
    dx[1] = x[2]
    dx[2] = -1/Q * x[2] - x[1] * Γ*sin(Ω*t + ϕ) + 24*V_0/k_c * σ^6 / (x[1] - d)^7 * (2 * σ^6/(x[1] - d)^6 - 1)
end


###TODO: add Lennard-Jones potential with Ingos damping term x^{-3}


"""
DMT model. Respects elastic properties via Hertzian contact model. 
"""
function f_DMT!(dx, x, p, t)
    Q, Ω, Γ, H, R, E, a_0, d, k_c, ϕ = p
    dx[1] = x[2]
    dx[2] = -1/Q * x[2] - x[1] + Γ*sin(Ω*t + ϕ)
    if x[1] < d - a_0
        dx[2] += H*R/(6*k_c * (d - x[1])^2)
    else
        dx[2] += H*R/(6 * a_0^2 * k_c) - 4/3*E*sqrt(R) * 1/(k_c) * (x[1] - (d - a_0))^(3/2) 
    end
end 


"""
DMT model with additional exponential decaying damping term, inspired by Haviland.
"""
function f_eDMT!(dx, x, p, t)
    γ_0, x_γ, ω_0, Q, Ω, Γ, H, R, E, a_0, d, k_c, ϕ = p
    dx[1] = x[2]
    dx[2] = -1/Q * x[2] - x[1] + Γ*sin(Ω*t + ϕ) + ω_0 * γ_0*exp((x[1] - d)/x_γ)*x[2]/k_c
    if x[1] < d - a_0
        dx[2] += H*R/(6*k_c * (d - x[1])^2)
    else
        dx[2] += H*R/(6 * a_0^2 * k_c) - 4/3*E*sqrt(R) * 1/(k_c) * (x[1] - (d - a_0))^(3/2) 
    end
end 


"""
DMR model with additional viscious damping term, only acting for x > a_0
"""
function f_vDMT!(dx, x, p, t)
    η, ω_0, Q, Ω, Γ, H, R, E, a_0, d, k_c, ϕ = p
    dx[1] = x[2]
    dx[2] = -1/Q * x[2] - x[1] + Γ*sin(Ω*t + ϕ) 
    if x[1] < d - a_0
        dx[2] += H*R/(6*k_c * (d - x[1])^2)
    else
        dx[2] += H*R/(6 * a_0^2 * k_c) - 4/3*E*sqrt(R) * 1/(k_c) * (x[1] - (d - a_0))^(3/2) + η/k_c * ω_0 * sqrt(R*(x[1] - a_0))*x[2]
    end
end 
