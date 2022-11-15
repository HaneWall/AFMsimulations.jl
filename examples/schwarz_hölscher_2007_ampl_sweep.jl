using AFMsimulations
using CairoMakie
CairoMakie.activate!(type="svg")

# parameters from schwarz and hölscher 2007
# force range?  

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

canti = Cantilever(R, Q, k, E_t, ν_t, f_0)
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

N = 50
f_high = 1.4e-9 
f_low = 1.240e-9 
f_fwd = collect(LinRange(f_low, f_high, N))

Δt = 0.005  # timestep Δτ in integration scheme, 2π/Δt steps per period (we use dimensionless time)

K = 10
Ω_low = 1.000
Ω_high = 1.003
Ω_fwd = collect(LinRange(Ω_low, Ω_high, K))


ampl_fwd = zeros(Float64, K, N)
ampl_bwd = zeros(Float64, K, N)

@inbounds for i in 1:K
    ampl_fwd[i, :] .= ampl_sweep(f_fwd, Ω_fwd[i], AFM, Δt)
    ampl_bwd[i, :] .= ampl_sweep(reverse(f_fwd), Ω_fwd[i], AFM, Δt)
end


fig = Figure(resolution=(800, 800), fontsize=24)
ax = Axis3(fig[1,1], xlabel=L"Ω", ylabel=L"F_{in}", zlabel=L"A_{out}")
scatter!(ax, Ω_fwd, f_fwd, ampl_fwd)
scatter!(ax, Ω_fwd, reverse(f_fwd), ampl_bwd)
fig


# # figure to represent the forward and backward sweep 
# fig = Figure(resolution=(800, 400), fontsize=24)
# ax = Axis(fig[1, 1], ylabel=L"A_{out}", xlabel=L"F_{in}")
# scatter!(ax,  f_fwd, ampl_fwd, marker=:utriangle, markersize=12, strokecolor=:black, strokewidth=1, color=(:black, 0.), label="Forward")
# scatter!(ax, reverse(f_fwd), ampl_bwd, marker=:dtriangle, markersize=12, strokecolor=:red, strokewidth=1, color=(:black, 0.), label="Backward")
# ylims!(ax, [5.5e-9, 10.5e-9])
# axislegend(ax, position=:cb)
# fig 
