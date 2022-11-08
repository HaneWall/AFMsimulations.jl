using AFMsimulations
using CairoMakie
using DifferentialEquations

CairoMakie.activate!(type = "svg")

# example sweep over the frequency for the DMT model 


canti = Cantilever(20.e-9, 300., 5., 1.e9, 0.3, 164973*2π)
probe = Sample(2.e-19, 1.e9, 0.3)
AFM = AFM_DMT_experiment(probe, canti, 0.3e-9, 8.5e-9)


force = 2.7e-9 # excitation force in Newton 
Δt = 0.005  # timestep Δτ in integration scheme, 2π/Δt steps per period (we use dimensionless time)
Ω = 1.003

T_sim = 4000
t_span= (0, T_sim)
u_0 = [0., 0.]
ϕ = 0.
Γ = force/canti.k
p = [AFM.tip.Q, Ω, Γ, AFM.sample.H, AFM.tip.R, AFM.eff_young, AFM.a_0, AFM.d, AFM.tip.k, ϕ]

prob = ODEProblem(f_DMT!, u_0, t_span, p)
sol = solve(prob, AutoTsit5(Rosenbrock23()), dt=Δt, adaptive=false) 

fig = Figure(resolution=(800, 800), fontsize=24)
ax = Axis(fig[1, 1])
T_period = 2π/Δt
timeslot_to_detect_ampl = ceil(Int, 2*T_period)
lines!(ax, sol[1, end-timeslot_to_detect_ampl:end], sol[2, end-timeslot_to_detect_ampl:end])
fig 
 