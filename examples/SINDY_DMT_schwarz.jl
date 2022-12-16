using AFMsimulations
using CairoMakie
using DifferentialEquations
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
ϕ = 0.




canti = Cantilever(R, Q, k, E_t, ν_t, f_0 * 2π)
probe = AFMsimulations.Sample(H, E_s, ν_s)
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
tspan = (0.0, 4800.0)
prob = ODEProblem(f_DMT!, u0, tspan, p)
sol = solve(prob, AutoTsit5(Rosenbrock23()), dt = 0.005, adaptive=false)

x_1_transient = sol[1, 1:400_000]
f_ts_transient = [get_tip_sample_interaction_dmt(x, p) for x in x_1_transient]

x_1 = sol[1, end-20_000:end]
f_ts = [get_tip_sample_interaction_dmt(x, p) for x in x_1] 

x_1_one = sol[1, end-1257:end]
f_ts_one = [get_tip_sample_interaction_dmt(x, p) for x in x_1_one]

fig = Figure(resolution = (1100, 1400))
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
using DataDrivenSparse

ts = sol.t[end-10_257:end]
X = sol[:, end-10_257:end] ./ d#.+ 0.01*randn(2, length(sol.t))
dx = 1/Δt * savitzky_golay_filter(X[1, :], 21, 3; deriv_order = 1, boundary_mode=:interpolation) 
dx2 = 1/Δt * savitzky_golay_filter(X[2, :], 21, 3; deriv_order = 1, boundary_mode=:interpolation)  
DX = transpose(hcat(dx, dx2))

fig = Figure(resolution=(400, 400))
ax = Axis(fig[2, 1])
lines!(ax, ts, X[1, :], linewidth = 2.)
lines!(ax, ts, X[2, :], linewidth = 2.)
lines!(ax, ts, dx)
lines!(ax, ts, dx2)
fig

@parameters t
@variables u(t)[1:2] c(t)[1]
@parameters w[1]

function forcing(u, p, t)
    return sin(t)
end

probdd = ContinuousDataDrivenProblem(X, ts, DX,
                                     U = (u, p, t) -> sin(t), 
                                     p=[1.])

function eff_function(l)
    IfElse.ifelse(l < d- a_0, inv((d - l)^2), (l - (d - a_0))^(3/2))
end

h = Num[polynomial_basis(u, 2); c]

basis = Basis(h, u, parameters = w, controls = c);
println(basis) # hide

sampler = DataProcessing(split = 0.8, shuffle = true, batchsize = 30)
#λs = exp10.(-60:0.1:-30)
opt = STLSQ(10e-9)
res = solve(probdd, basis, opt, options = DataDrivenCommonOptions(data_processing = sampler, digits = 2))

system = get_basis(res)
params = get_parameter_map(system)
println(system) # hide
println(params) # hide


t = Symbolics.unwrap(get_iv(system))
subs_control = Dict(
    c[1] => forcing(u, p, t) 
)


eqs = map(equations(system)) do eq
    eq.lhs ~ substitute(eq.rhs, subs_control)
end

@named sys = ODESystem(
    eqs,
    get_iv(system),
    states(system),
    parameters(system)
    )

x0 = [u[1] => X[1, 1], u[2] => X[2, 1]]
ps = get_parameter_map(system)

t_span = (ts[1], ts[end])
ode_prob = ODEProblem(sys, x0, t_span, ps)
estimate = solve(ode_prob, Tsit5(), dt=Δt, adaptive=false)


fig = Figure(resolution=(800, 800), fontsize = 24, font=("CMU Serif", ))
ax = Axis(fig[1, 1], xlabel=L"t")
lines!(ax, ts, X[1, :], label=L"x")
lines!(ax, ts, X[2, :], label=L"\dot{x}")
lines!(ax, estimate.t, estimate[1, :], linestyle=:dash, label=L"x_{Sindy}")
lines!(ax, estimate.t, estimate[2, :], linestyle=:dash, label=L"\dot{x}_{Sindy}")
lines!(ax, ts, 1e-8 .* sin.(ts))
axislegend()
fig

# polys = polynomial_basis(u, 2)
# push!(polys, sin.(u[1]))
# push!(polys, cos.(u[1]))
# push!(polys, sin.(u[1])^2)
# push!(polys, cos.(u[1])^2)
# push!(polys, sin.(u[1]) .* u[3:4]...)
# push!(polys, sin.(u[1]) .* u[3:4] .^ 2...)
# push!(polys, sin.(u[1]) .* cos.(u[1])...)
# push!(polys, sin.(u[1]) .* cos.(u[1]) .* u[3:4]...)
# push!(polys, sin.(u[1]) .* cos.(u[1]) .* u[3:4] .^ 2...)

# h = Num[ ;c]


# t_span = (4993.715000085791, 5000.)

# basis = Basis(implicits, u, implicits = du, controls = x, iv = t);
