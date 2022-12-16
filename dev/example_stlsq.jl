using DataDrivenDiffEq
using DataDrivenSparse
using LinearAlgebra
using OrdinaryDiffEq
using CairoMakie
CairoMakie.activate!(type="svg")
using AFMsimulations

function pendulum(u, p, t)
    x = u[2]
    y = -9.81sin(u[1]) - 0.3u[2]^3 - 3.0 * cos(u[1]) - 10.0 * exp(-((t - 5.0) / 5.0)^2)
    return [x; y]
end

u0 = [0.99π; -1.0]
tspan = (0.0, 30.0)
Δt = 0.01
prob = ODEProblem(pendulum, u0, tspan)
sol = solve(prob, Tsit5(), dt = Δt, adaptive = false)


X = sol[:, :] #.+ 0.01*randn(2, length(sol.t))
dx = 1/Δt * savitzky_golay_filter(X[1, :], 71, 3; deriv_order = 1, boundary_mode=:nearest) 
dx2 = 1/Δt * savitzky_golay_filter(X[2, :], 71, 3; deriv_order = 1, boundary_mode=:nearest)  
DX = transpose(hcat(dx, dx2))

#using CairoMakie
fig = Figure(resolution=(400, 800))
ax = Axis(fig[1, 1])
lines!(ax, sol.t, sol[1, :])
lines!(ax, sol.t, sol[2, :])
ax2 = Axis(fig[2, 1])
lines!(ax2, sol.t, X[1, :], linewidth = 2.)
lines!(ax2, sol.t, X[2, :], linewidth = 2.)
lines!(ax2, sol.t, dx)
lines!(ax2, sol.t, dx2)
fig


ts = sol.t;

function forcing(u, p, t)
    return exp(-((t - 5.0) / 5.0)^2)
end

prob = ContinuousDataDrivenProblem(X, ts, DX,
                                   U = (u, p, t) -> forcing(u, p, t),
                                   p = ones(2))
@parameters t
@variables u(t)[1:2] c(t)[1:1]
@parameters w[1:2]
u = collect(u)
c = collect(c)
w = collect(w)

h = Num[sin.(w[1] .* u[1]); cos.(w[2] .* u[1]); polynomial_basis(u, 4); c]

basis = Basis(h, u, parameters = w, controls = c);
println(basis) # hide

sampler = DataProcessing(split = 0.8, shuffle = true, batchsize = 30)
λs = exp10.(-10:0.1:0)
opt = STLSQ(0.2)
res = solve(prob, basis, opt, options = DataDrivenCommonOptions(data_processing = sampler, digits = 2, selector = bic))

system = get_basis(res)
params = get_parameter_map(system)
println(system) # hide
println(params) # hide


t = Symbolics.unwrap(get_iv(system))
subs_control = Dict(
    c[1] => exp(-((t - 5.0) / 5.0)^2)
)


eqs = map(equations(system)) do eq
    eq.lhs ~ substitute(eq.rhs, subs_control)
end

@named sys = ODESystem(
    eqs,
    get_iv(system),
    states(system),
    parameters(system)
    );

x0 = [u[1] => u0[1], u[2] => u0[2]]
ps = get_parameter_map(system)

ode_prob = ODEProblem(sys, x0, tspan, ps)
estimate = solve(ode_prob, Tsit5(), dt=Δt, adaptive=false)


fig = Figure(resolution=(800, 800), fontsize = 24, font=("CMU Serif", ))
ax = Axis(fig[1, 1], xlabel=L"t")
lines!(ax, sol.t, X[1, :], label=L"x")
lines!(ax, sol.t, X[2, :], label=L"\dot{x}")
lines!(ax, estimate.t, estimate[1, :], linestyle=:dash, label=L"x_{Sindy}")
lines!(ax, estimate.t, estimate[2, :], linestyle=:dash, label=L"\dot{x}_{Sindy}")
axislegend()
fig
