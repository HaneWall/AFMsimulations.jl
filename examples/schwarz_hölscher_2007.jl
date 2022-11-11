using AFMsimulations
using CairoMakie
CairoMakie.activate!(type = "svg")

# parameters from schwarz and hölscher 2007

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


canti = Cantilever(R, Q, k, E_t, ν_t, f_0 * 2π)
probe = Sample(H, E_s, ν_s)
AFM = AFM_DMT_experiment(probe, canti, a_0, d)

# # Plot force curve 
# x = collect(LinRange(-5.e-9, 10.e-9, 500))
# F = [force_distance(x[i], AFM) for i in eachindex(x)]
# fig = Figure(resolution=(800, 800))
# axf = Axis(fig[1, 1], ylabel=L"F[N]", xlabel=L"x[m]")
# lines!(axf, x, F)
# #lines!(axf, x, -k * x)
# fig


N = 170 # number of omegas for one forward or backward sweep  
Ω_low = 0.9975
Ω_high = 1.0025
Ω_fwd = collect(LinRange(Ω_low, Ω_high, N))
Δt = 0.005  # timestep Δτ in integration scheme, 2π/Δt steps per period (we use dimensionless time)


ampl_fwd = freq_sweep(force, Ω_fwd, AFM, Δt)
ampl_bwd = freq_sweep(force, reverse(Ω_fwd), AFM, Δt)

# figure to represent the forward and backward sweep 
fig = Figure(resolution=(800, 400), fontsize=24)
ax = Axis(fig[1, 1], ylabel=L"A", xlabel=L"f[Hz]")
scatter!(ax, Ω_fwd .* 300000, ampl_fwd, marker=:utriangle, markersize=12, strokecolor=:black, strokewidth=1, color=(:black, 0.), label="Forward")
scatter!(ax, reverse(Ω_fwd) .* 300000 , ampl_bwd, marker=:dtriangle, markersize=12, strokecolor=:red, strokewidth=1, color=(:black, 0.), label="Backward")
xlims!(ax, [Ω_low * 300000, Ω_high * 300000])
ylims!(ax, [5.5e-9, 10.5e-9])
axislegend(ax, position=:cb)
fig 

