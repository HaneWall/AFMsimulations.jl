using AFMsimulations
using CairoMakie
CairoMakie.activate!(type="svg")

# Let us start by getting a normal uncontrolled frequency sweep for a vLJ-model

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
barke_probe = Sample(H, E_s, ν_s)
AFM = AFM_vLJ_experiment(barke_probe, barke_canti, σ, δx, V_0, γ, d)

N = 180         # number of omegas per forward/backward sweep
Ω_low = 0.965   # big Omega is the ratio ω_drive/ω_0
Ω_high = 1.035
Ω_fwd = collect(LinRange(Ω_low, Ω_high, N))
Δt = 0.005      # timestep Δτ in integration scheme, 2π/Δt steps per period (we use dimensionless time)

# uncontrolled data 
ampl_fwd, phi_fwd = freq_sweep_stat(force, Ω_fwd, AFM, Δt)
ampl_bwd, phi_bwd = freq_sweep_stat(force, reverse(Ω_fwd), AFM, Δt)

# figure to represent the uncontrolled forward and backward sweep 
fig = Figure(resolution=(800, 1300), fontsize=24, fonts = (; regular = "CMU Serif"))
ax = Axis(fig[1, 1], ylabel=L"A", xlabel=L"f[Hz]")
scatter!(ax, Ω_fwd .* f_0, ampl_fwd, marker=:utriangle, markersize=12, strokecolor=:black, strokewidth=1, color=(:black, 0.), label="Forward")
scatter!(ax, reverse(Ω_fwd) .* f_0 , ampl_bwd, marker=:dtriangle, markersize=12, strokecolor=:red, strokewidth=1, color=(:black, 0.), label="Backward")
xlims!(ax, [Ω_low * f_0, Ω_high * f_0])
axislegend(ax, position=:cb)
ax2 = Axis(fig[2, 1], ylabel=L"\varphi", xlabel=L"f[Hz]") 
scatterlines!(ax2, Ω_fwd .* f_0, phi_fwd, marker=:utriangle, markersize=12, strokecolor=:black, strokewidth=1, color=(:black, 0.), label="Forward")
scatterlines!(ax2, reverse(Ω_fwd) .* f_0 , phi_bwd, marker=:dtriangle, markersize=12, strokecolor=:red, strokewidth=1, color=(:red, 0.), label="Backward")
xlims!(ax2, [Ω_low * f_0, Ω_high * f_0])
axislegend(ax2, position=:lt)
fig 

# anna PseudoCode

φ_target = 
φ

