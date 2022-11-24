using StaticArrays, BenchmarkTools, DifferentialEquations, AFMsimulations

function f_vLJ_static(u, p, t)
    σ, δx, V_0, γ, ω_0, Q, Ω, Γ, d, k_c, ϕ = p
    # define some helping variables (different coordinate system)
    h_x = (u[1] - d)^2 + (δx)^2
    dx = u[2]
    dy = -1/Q * u[2] - u[1] + Γ*sin(Ω*t + ϕ) - 12*V_0/(k_c * sqrt(h_x)) * ((σ^2 / h_x)^6 - (σ^2 / h_x)^3)  - u[2] * ω_0 * γ/(d - u[1])^3
    SA[dx, dy]
end


function f_vLJ!(dx, x, p, t)
    σ, δx, V_0, γ, ω_0, Q, Ω, Γ, d, k_c, ϕ = p
    # define some helping variables (different coordinate system)
    h_x = (x[1] - d)^2 + (δx)^2
    dx[1] = x[2]
    dx[2] = -1/Q * x[2] - x[1] + Γ*sin(Ω*t + ϕ) - 12*V_0/(k_c * sqrt(h_x)) * ((σ^2 / h_x)^6 - (σ^2 / h_x)^3)  - x[2] * ω_0 * γ/(d - x[1])^3
end


### Barke parameters

d = 100.e-9             # distance equilibrium cantilever to surface
Q = 400.                # quality factor of the spring 
k = 0.7                 # spring constant 
f_0 = 50.e3             # eigenfrequency cantilever
A_free = 250.e-9        # free amplitude
force = k * A_free/Q

V_0 = 4.2e-18           # Epsilon parameter bei Ingo, ACHTUNG er nutzt nur nm Skala --> 4.2 zu 4.2e-18
σ = 2.8e-9              # Ingo (zz)
δx = 0.5e-9             # softening parameter
γ = 1.e-31              # 1.e-4 in Nm nm^2 s


### some usual parameters for a simulation 
p = [σ, δx, V_0, γ, 2π*f_0, Q, 1., force/k, d, k, 0.]

### define staticarray problem
u0 = SA[0.0, 0.0]
tspan = (0., 8000.)
stat_prob = ODEProblem(f_vLJ_static, u0, tspan, p)
@benchmark solve(stat_prob, AutoTsit5(Rosenbrock23()), dt = 0.005, adaptive=false)


### define in-place problem 
u0 = [0., 0.]
tspan = (0., 8000.)
ip_prob = ODEProblem(f_vLJ!, u0, tspan, p)
@benchmark solve(ip_prob, AutoTsit5(Rosenbrock23()), dt = 0.005, adaptive=false)

"""
Static Arrays always outperform in-place evaluations for small ODE's. Speedup: 1x..3x
"""