using DifferentialEquations, AFMsimulations, BenchmarkTools, StaticArrays

"""
In order to get more convenient and faster computations it would be quite nice,
to use the integrator interface from DifferentialEquations.jl. This would allow us
to use Callbacks for the integrations. 
"""

function f_vLJ_static(u, p, t)
    σ, δx, V_0, γ, ω_0, Q, Ω, Γ, d, k_c, ϕ = p
    # define some helping variables (different coordinate system)
    h_x = (u[1] - d)^2 + (δx)^2
    dx = u[2]
    dy = -1/Q * u[2] - u[1] + Γ*sin(Ω*t + ϕ) - 12*V_0/(k_c * sqrt(h_x)) * ((σ^2 / h_x)^6 - (σ^2 / h_x)^3)  - u[2] * ω_0 * γ/(d - u[1])^3
    SA[dx, dy]
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


u0 = SA[0.0, 0.0]
tspan = (0., 8000.)
stat_prob = ODEProblem(f_vLJ_static, u0, tspan, p)


### define staticarray problem (with solve setup - no callback if periodic orbit is reached)
@benchmark solve(stat_prob, AutoTsit5(Rosenbrock23()), dt = 0.005, adaptive=false)


### define integrator setup 
@btime integrator = init(stat_prob, AutoTsit5(Rosenbrock23()), dt = 0.005, adaptive=false)
@benchmark solve!(integrator)

"""
INSANE speedup, we are in in NANOSECONDS-regime now. 
"""