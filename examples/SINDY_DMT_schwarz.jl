using AFMsimulations
using CairoMakie
using DifferentialEquations
#using DataDrivenDiffEq 
#using ModelingToolkit
using LinearAlgebra
CairoMakie.activate!(type="svg")


force = 1.333e-9
H = 2.e-19
a_0 = 0.3e-9
d = 8.5e-9

E_s = 1.e9
ν_s = 0.3 

E_t = 130.e9
ν_t = 0.3
Q = 300.
k = 40.
f_0 = 300000.
R = 10.e-9
Ω = 1. 
Γ = force / k 




canti = Cantilever(R, Q, k, E_t, ν_t, f_0 * 2π)
probe = AFMSimulations.Sample(H, E_s, ν_s)
AFM = AFM_DMT_experiment(probe, canti, a_0, d)

E = eff_young_module(probe, canti)

p = Q, Ω, Γ, H, R, E, a_0, d, k, ϕ 


function get_tip_sample_interaction_dmt(u1, p)
    Q, Ω, Γ, H, R, E, a_0, d, k, ϕ  = p
    if u1 < d - a_0
        f_ts = H*R/(6*k * (d - u1)^2)
    else
        f_ts = H*R/(6 * a_0^2 * k) - 4/3*E*sqrt(R) * 1/(k) * (u1 - (d - a_0))^(3/2) 
    end
    return f_ts 
end


u0 = [0.; 0.]
tspan = (0.0, 5000.0)
prob = ODEProblem(f_DMT!, u0, tspan, p)
sol = solve(prob, AutoTsit5(Rosenbrock23()), dt = 0.005, adaptive=false)

x_1_transient = sol[1, 1:400_000]
f_ts_transient = [get_tip_sample_interaction_dmt(x, p) for x in x_1_transient]

x_1 = sol[1, end-20_000:end]
f_ts = [get_tip_sample_interaction_dmt(x, p) for x in x_1] 

x_1_one = sol[1, end-1257:end]
f_ts_one = [get_tip_sample_interaction_dmt(x, p) for x in x_1_one]

fig = Figure(resolution = (1100, 1200))
ax = Axis3(fig[1, 1])
lines!(ax, sol[1, 1:400_000], sol[2, 1:400_000], sol.t[1:400_000], linewidth = .3)
ax_ts = Axis(fig[1, 2])
lines!(ax_ts, sol.t[1:400_000], f_ts_transient)
ax_stab = Axis3(fig[2, 1])
lines!(ax_stab, sol[1, end-20_000:end], sol[2, end-20_000:end], sol.t[end-20_000:end]; color=sol.t[end-20_000:end], colormap=:Spectral_11, linewidth = 1.)
ax_stab_ts = Axis(fig[2, 2])
lines!(ax_stab_ts, sol.t[end-20_000:end], f_ts; color=sol.t[end-20_000:end], colormap=:Spectral_11)
ax_one = Axis(fig[3, 1])
lines!(ax_one, sol[1, end-1257:end], sol[2, end-1257:end];color=sol.t[end-1257:end], colormap=:hsv, linewidth = 2.)
ax_stab_ts_one = Axis(fig[3, 2])
lines!(ax_stab_ts_one, sol.t[end-1257:end], f_ts_one; color=sol.t[end-1257:end], colormap=:hsv)
fig


using DataDrivenDiffEq 
using ModelingToolkit

@parameters t
@variables u[1:2] c[1:1]

polys = polynomial_basis(u, 2)
push!(polys, sin.(u[1]))
push!(polys, cos.(u[1]))
push!(polys, sin.(u[1])^2)
push!(polys, cos.(u[1])^2)
push!(polys, sin.(u[1]) .* u[3:4]...)
push!(polys, sin.(u[1]) .* u[3:4] .^ 2...)
push!(polys, sin.(u[1]) .* cos.(u[1])...)
push!(polys, sin.(u[1]) .* cos.(u[1]) .* u[3:4]...)
push!(polys, sin.(u[1]) .* cos.(u[1]) .* u[3:4] .^ 2...)

h = Num[ ;c]


t_span = (4993.715000085791, 5000.)

basis = Basis(implicits, u, implicits = du, controls = x, iv = t);
