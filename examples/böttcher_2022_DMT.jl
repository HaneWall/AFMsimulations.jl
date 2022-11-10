using AFMsimulations
using CairoMakie
CairoMakie.activate!(type = "svg")

# try to recreate sweeps from Lukas Böttcher (ppt)


H = 2.96e-19    # Si-Si interaction Hamaker
a_0 = 0.3e-9
d = 24.e-9

E_s = 1.e9
ν_s = 0.3 

E_t = 130.e9    # Si Young's module 
ν_t = 0.3
Q = 300.
k = 7.
f_0 = 54000.
R = 10.e-9

A_free = 25.e-9
force = k * A_free/Q


canti = Cantilever(R, Q, k, E_t, ν_t, f_0 * 2π)
probe = Sample(H, E_s, ν_s)
AFM = AFM_DMT_experiment(probe, canti, a_0, d)

N = 150
Ω_low = 0.975 
Ω_high = 1.025 
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
#ylims!(ax, [5.5e-9, 10.5e-9])
axislegend(ax, position=:cb)
fig 


