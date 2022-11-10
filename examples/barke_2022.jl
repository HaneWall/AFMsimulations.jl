using AFMsimulations
using CairoMakie
CairoMakie.activate!(type="svg")
# Trying to reproduce Ingo Sweep, powerpoint page 20 (LennardJones.pptx - Unibox)

d = 100.e-9
Q = 400. 
k = 0.7
f_0 = 50.e3
#γ = 2.5e-14 # damping term for x^{-3} dependence 
A_free = 250.e-9
force = k * A_free/Q

H = 2.9e-19
V_0 = 4.2e-18
σ = 2.8e-9
δx = 0.5e-9
γ = 0. #-1.e-4




#some parameters that dont matter in LJ-sweep anyway (need for construction of objects)
R = 10.e-9 
E_t = 130.e9
ν_t = 0.3

E_s = 1.e-9
ν_s = 0.3

barke_canti = Cantilever(R, Q, k, E_t, ν_t, f_0)
barke_probe = Sample(H, E_s, ν_s)
AFM = AFM_vLJ_experiment(barke_probe, barke_canti, σ, δx, V_0, γ, d)

# #let us plot the force distance relation (also seen in Ingos Figure)
# x = collect(LinRange(2.e-9, 6.e-9, 200))
# F = [force_distance(x[i], AFM) for i in eachindex(x)]
# fig = Figure(resolution=(800, 800))
# axf = Axis(fig[1, 1], ylabel=L"F[N]", xlabel=L"x[m]")
# lines!(axf, x, F)
# ylims!(axf, [-11.e-9, 31.e-9])
# fig

# # map z-d --> x , z = x + d

# x_h = collect(LinRange(-5.e-9, 6.e-9, 200))
# d_h = 8.e-9# 10.e-9 
# δx = 0.5e-9

# diffe(x) = x^2 + δx^2
# vdW(x) = (σ^2/diffe(x))^3
# paul(x) = (σ^2/diffe(x))^6
# f(x ) = 12*V_0/sqrt(diffe(x)) * (paul(x) - vdW(x))
# h_f = [f(i - d_h) for i in x_h]
# fig = Figure(resolution=(800, 800))
# axf = Axis(fig[1, 1], ylabel=L"F[N]", xlabel=L"x[m]")
# lines!(axf, x_h, h_f)
# ylims!(axf, [-11e-9, 30e-9])
# fig



N = 75
Ω_low = 0.965 
Ω_high = 1.035 
Ω_fwd = collect(LinRange(Ω_low, Ω_high, N))
Δt = 0.005  # timestep Δτ in integration scheme, 2π/Δt steps per period (we use dimensionless time)


ampl_fwd = freq_sweep(force, Ω_fwd, AFM, Δt)
ampl_bwd = freq_sweep(force, reverse(Ω_fwd), AFM, Δt)

# figure to represent the forward and backward sweep 
fig = Figure(resolution=(800, 400), fontsize=24)
ax = Axis(fig[1, 1], ylabel=L"A", xlabel=L"f[Hz]")
scatter!(ax, Ω_fwd .* f_0, ampl_fwd, marker=:utriangle, markersize=12, strokecolor=:black, strokewidth=1, color=(:black, 0.), label="Forward")
scatter!(ax, reverse(Ω_fwd) .* f_0 , ampl_bwd, marker=:dtriangle, markersize=12, strokecolor=:red, strokewidth=1, color=(:black, 0.), label="Backward")
xlims!(ax, [Ω_low * f_0, Ω_high * f_0])
#ylims!(ax, [50e-9, 10.5e-9])
axislegend(ax, position=:cb)
fig 