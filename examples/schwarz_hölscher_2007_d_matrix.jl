using AFMsimulations
using CairoMakie
using Printf
CairoMakie.activate!(type = "svg")

# parameters from schwarz and hölscher 2007, small datamatrix over d 

force = 1.333e-9
H = 2.e-19
a_0 = 0.3e-9
d = [8.0e-9, 8.5e-9, 9.e-9, 9.5e-9]

E_s = 1.e9
ν_s = 0.3 

E_t = 130.e9
ν_t = 0.3
Q = 300.
k = 40.
f_0 = 300000.
R = 10.e-9


canti = Cantilever(R, Q, k, E_t, ν_t, f_0)
probe = Sample(H, E_s, ν_s)
AFM_exps = [AFM_DMT_experiment(probe, canti, a_0, d[i]) for i in 1:4]

N = 130 # number of omegas for one forward or backward sweep  
Ω_low = 0.997
Ω_high = 1.003
Ω_fwd = collect(LinRange(Ω_low, Ω_high, N))
Δt = 0.005  # timestep Δτ in integration scheme, 2π/Δt steps per period (we use dimensionless time)
ampl_fwd = zeros(Float64, 4, N)
ampl_bwd = zeros(Float64, 4, N)

for i in 1:4
    @inbounds ampl_fwd[i, :] .= freq_sweep(force, Ω_fwd, AFM_exps[i], Δt) 
    @inbounds ampl_bwd[i, :] .= freq_sweep(force, reverse(Ω_fwd), AFM_exps[i], Δt)
end


# figure to represent the forward and backward sweeps
fig = Figure(resolution=(1200, 800), fontsize=24)

ax1 = Axis(fig[1, 1], ylabel=L"A[nm]", xlabel=L"f_{Drive}[kHz]")
scatter!(ax1, Ω_fwd .* 300, ampl_fwd[1, :] .*1.e9, marker=:utriangle, markersize=12, strokecolor=:black, strokewidth=1, color=(:black, 0.), label="Forward")
scatter!(ax1, reverse(Ω_fwd) .* 300 , ampl_bwd[1, :].*1.e9, marker=:dtriangle, markersize=12, strokecolor=:red, strokewidth=1, color=(:black, 0.), label="Backward")
xlims!(ax1, [299., 301.])
ylims!(ax1, [5.5, 10.5])
text!(
    "d=" .* string.(d[1]*1.e9) .*"nm",
    position = Point2f(299.5, 10.),
    align = (:right, :baseline),
)
axislegend(ax1, position=:cb)


ax2 = Axis(fig[1, 2], ylabel=L"A[nm]", xlabel=L"f_{Drive}[kHz]")
scatter!(ax2, Ω_fwd .* 300, ampl_fwd[2, :].*1.e9, marker=:utriangle, markersize=12, strokecolor=:black, strokewidth=1, color=(:black, 0.), label="Forward")
scatter!(ax2, reverse(Ω_fwd) .* 300 , ampl_bwd[2, :].*1.e9, marker=:dtriangle, markersize=12, strokecolor=:red, strokewidth=1, color=(:black, 0.), label="Backward")
xlims!(ax2, [299., 301.])
ylims!(ax2, [5.5, 10.5])
text!(
    "d=" .* string.(d[2]*1.e9) .*"nm",
    position = Point2f(299.5, 10.),
    align = (:right, :baseline),
)
axislegend(ax2, position=:cb)

ax3 = Axis(fig[2, 1], ylabel=L"A[nm]", xlabel=L"f_{Drive}[kHz]")
scatter!(ax3, Ω_fwd .* 300, ampl_fwd[3, :] .* 1.e9, marker=:utriangle, markersize=12, strokecolor=:black, strokewidth=1, color=(:black, 0.), label="Forward")
scatter!(ax3, reverse(Ω_fwd) .* 300 , ampl_bwd[3, :].*1.e9, marker=:dtriangle, markersize=12, strokecolor=:red, strokewidth=1, color=(:black, 0.), label="Backward")
xlims!(ax3, [299., 301.])
ylims!(ax3, [5.5, 10.5])
text!(
    "d=" .* string.(d[3]*1.e9) .*"nm",
    position = Point2f(299.5, 10.),
    align = (:right, :baseline),
)
axislegend(ax3, position=:cb)

ax4 = Axis(fig[2, 2], ylabel=L"A[nm]", xlabel=L"f_{Drive}[kHz]")
scatter!(ax4, Ω_fwd .* 300, ampl_fwd[4, :] .*1.e9, marker=:utriangle, markersize=12, strokecolor=:black, strokewidth=1, color=(:black, 0.), label="Forward")
scatter!(ax4, reverse(Ω_fwd) .* 300 , ampl_bwd[4, :] .* 1.e9, marker=:dtriangle, markersize=12, strokecolor=:red, strokewidth=1, color=(:black, 0.), label="Backward")
xlims!(ax4, [299., 301.])
ylims!(ax4, [5.5, 10.5])
text!(
    "d=" .* string.(d[4]*1.e9) .*"nm",
    position = Point2f(299.5, 10.),
    align = (:right, :baseline),
)
axislegend(ax4, position=:cb)

fig 

