# import the package we just wrote in src 
using AFMsimulations
using CairoMakie
CairoMakie.activate!(type = "svg")

# example sweep over the frequency for the DMT model 


canti = Cantilever(20.e-9, 300., 5., 109.e9, 0.3, 164973*2π)
probe = Sample(2.e-19, 1.e9, 0.3)
AFM = AFM_DMT_experiment(probe, canti, 0.3e-9, 8.5e-9)


N = 100 # number of omegas for one forward or backward sweep  
Ω_low = 0.98
Ω_high = 1.03
Ω_fwd = collect(LinRange(Ω_low, Ω_high, N))
force = 0.15e-9 # excitation force in Newton 
Δt = 0.005  # timestep Δτ in integration scheme, 2π/Δt steps per period (we use dimensionless time)


ampl_fwd = freq_sweep(force, Ω_fwd, AFM, Δt)
ampl_bwd = freq_sweep(force, reverse(Ω_fwd), AFM, Δt)

# figure to represent the forward and backward sweep 
fig = Figure(resolution=(800, 800), fontsize=24)
ax = Axis(fig[1, 1], ylabel=L"A", xlabel=L"ω_{Drive}/ω_0")

scatter!(ax, Ω_fwd, ampl_fwd, marker=:circle, markersize=12, strokecolor=:black, strokewidth=1, color=(:black, 0.), label="Forward")
scatter!(ax, reverse(Ω_fwd), ampl_bwd, marker=:circle, markersize=12, strokecolor=:blue, strokewidth=1, color=(:black, 0.), label="Backward")
xlims!(ax, [Ω_low, Ω_high])
axislegend(ax, position=:cb)
fig 