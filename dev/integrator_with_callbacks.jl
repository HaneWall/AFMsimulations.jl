using DifferentialEquations, StaticArrays, AFMsimulations, BenchmarkTools

"""
Idea: Terminate integration if the solution is periodic. That is x_vec(t + T) = x_vec(t). 
Implementation: create Matrix: 
A = [   
    x_1(t-T) x_2(t-T);
    x_1(t)   x_2(t)
]

--> This describes a saving callback. 
--> push! :new array at bottom x_vec(t)
--> popfirst! :delete Datapoint from before x_vec(t-2T)
--> norm(A[1, :] - A[2, :]) < ϵ --> Converged onto periodic orbit. This is like a Poincare-criterium. 
--> terminate the integrator
--> read ϕ, and the amplitude --> save both values
--> create new Integrator, Even better: affect the integrator 
"""

function Poincare_pass(integrator, abstol, T_period, t_min)
    popfirst!(M)
    push!(M, integrator.u)
    if norm(M[1, :] - M[2, :]) < abstol && integrator.t >= t_min
        return true
    else
        return false 
    end
end

function 

function TerminatePeriodicOrbit(abstol = 1e-10, T_period = 2π; min_t = nothing)
    condition = (u, t, integrator) -> Poincare_pass(integrator, abstol, T_period, min_t)
    affect! = (integrator) -> terminate!(integrator)
    return DiscreteCallback(condition, affect!; save_positions = (true, false))
end


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
Ω = 1. 
p = [σ, δx, V_0, γ, 2π*f_0, Q, Ω, force/k, d, k, 0.]


Δt = 0.005
u0 = SA[0.0, 0.0]
tspan = (0., 8000.)
stat_prob = ODEProblem(f_vLJ_static, u0, tspan, p)
T = Ω * 2π

integrator = init(stat_prob, AutoTsit5(Rosenbrock23()), dt=Δt, adaptive=false, cb=TerminatePeriodicOrbit(abstol = 1e-10, T_period = 2π))

solve!(integrator)