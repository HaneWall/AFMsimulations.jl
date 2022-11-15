"""
Lennard-Jones model, smooth interaction potential.
"""
function f_LJ!(dx, x, p, t)
    σ, V_0, Q, Ω, Γ, d, k_c, ϕ = p
    dx[1] = x[2]
    dx[2] = -1/Q * x[2] - x[1] + Γ*sin(Ω*t + ϕ) + 24*V_0/k_c * σ^6 / (x[1] - d)^7 * (2 * σ^6/(x[1] - d)^6 - 1)
end


"""
Lennard-Jones model, smooth interaction, soft sphere model to avoid
singularities at the surface of the sample. Introduces damping via (x - d)^{-3}
term. See I. Barke implementation. 
"""
function f_vLJ!(dx, x, p, t)
    σ, δx, V_0, γ, ω_0, Q, Ω, Γ, d, k_c, ϕ = p
    # define some helping variables (different coordinate system)
    h_x = (x[1] - d)^2 + (δx)^2
    dx[1] = x[2]
    dx[2] = -1/Q * x[2] - x[1] + Γ*sin(Ω*t + ϕ) - 12*V_0/(k_c * sqrt(h_x)) * ((σ^2 / h_x)^6 - (σ^2 / h_x)^3)  - x[2] * ω_0 * γ/(d - x[1])^3
end



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
Autonomous DMT model. Respects elastic properties via Hertzian contact model. 
"""
function f_DMT_auto!(dx, x, p, t)
    Q, Ω, Γ, H, R, E, a_0, d, k_c, ϕ = p
    dx[1] = x[2]
    dx[2] = -1/Q * x[2] - x[1] + Γ*sin(Ω*x[3] + ϕ)
    dx[3] = (1.)%2π 
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
    dx[2] = -1/Q * x[2] - x[1] + Γ*sin(Ω*t + ϕ) - ω_0 * γ_0*exp((x[1] - d)/x_γ)*x[2]/k_c
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
