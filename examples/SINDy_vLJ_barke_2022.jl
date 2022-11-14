using AFMsimulations
using CairoMakie
using DataDrivenDiffEq 
using ModelingToolkit

using LinearAlgebra
CairoMakie.activate!(type="svg")
# Trying to reproduce Ingo Sweep, powerpoint page 22 (LennardJones.pptx - Unibox)

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


# some parameters that dont matter in LJ-sweep anyway (need for construction of objects)
# dummy parameters in this simulation
H = 2.9e-19
R = 10.e-9 
E_t = 130.e9
ν_t = 0.3
E_s = 1.e-9
ν_s = 0.3

# construction of the simulation 
barke_canti = Cantilever(R, Q, k, E_t, ν_t, f_0)
barke_probe = AFMsimulations.Sample(H, E_s, ν_s)
AFM = AFM_vLJ_experiment(barke_probe, barke_canti, σ, δx, V_0, γ, d)

# work on the reosnance frequency
Ω = 1. 
ϕ = 0.
Γ = force/k
p = [AFM.σ, AFM.δx, AFM.V_0, AFM.γ, 2π*AFM.tip.f_0, AFM.tip.Q, Ω, Γ, AFM.d, AFM.tip.k, ϕ]

u0 = [0.; 0.]
tspan = (0.0, 5000.0)
prob = ODEProblem(f_vLJ!, u0, tspan, p)
sol = solve(prob, AutoTsit5(Rosenbrock23()), dt = 0.005, adaptive=false)

fig = Figure(resolution=(600, 600))
ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"∂x")
lines!(ax, sol[1, :], sol[2, :])
fig 


## SINDy reconstruction without control 
# Create Datamatrix X
X = sol[:,:] #+ 0.2 .* randn(size(sol)) 
ts = sol.t


#### Start the automatic discovery
ddprob = ContinuousDataDrivenProblem(sol)
@variables u[1:2] c[1:1]
@parameters w[1:2]
u = collect(u)
c = collect(c)
w = collect(w)

h = Num[sin.(w[1].*u[1]);cos.(w[2].*u[1]); polynomial_basis(u, 5); c]

basis = Basis(h, u, parameters = w, controls = c);



@variables t x(t) y(t) 
u = [x;y;z]
basis = Basis(polynomial_basis(u, 5), u, iv = t)
opt = STLSQ(exp10.(-5:0.1:-1))
ddsol = solve(ddprob, basis, opt, normalize = true)
print(ddsol, Val{true})
