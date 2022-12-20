using DifferentialEquations, DiffEqCallbacks, StaticArrays, Parameters, BenchmarkTools, LinearAlgebra, Statistics, CairoMakie




function f_vLJ_static(u, p, t)
    @unpack σ, δx, V_0, γ, ω_0, Q, Ω, force, d, k, ϕ = p
    # define some helping variables (different coordinate system)
    h_x = (u[1] - d)^2 + (δx)^2
    dx = u[2]
    dy = -1/Q * u[2] - u[1] + force*sin(Ω*t + ϕ) - 12*V_0/(k * sqrt(h_x)) * ((σ^2 / h_x)^6 - (σ^2 / h_x)^3)  - u[2] * ω_0 * γ/(d - u[1])^3
    SA[dx, dy]
end


"""
DMT model with StaticArrays.
"""
function f_DMT(x, p, t)
    @unpack Q, Ω, Γ, H, R, E, a_0, d, k, ϕ = p
    dx = x[2]
    dy = -1/Q * x[2] - x[1] + Γ*sin(Ω*t + ϕ)
    if x[1] < d - a_0
        dy += H*R/(6*k * (d - x[1])^2)
    else
        dy += H*R/(6 * a_0^2 * k) - 4/3*E*sqrt(R) * 1/(k) * (x[1] - (d - a_0))^(3/2) 
    end
    SA[dx, dy]
end 

function lsq_regression(X, y, bases)
    B = [b(x) for x in X, b in bases]
    θ = (B'B)\B'y
    return θ
end

function sinusoidal_bases_1d(j, k, Ω, a, b) 
    T = b[j] - a[j]
    bases = Function[x->1/2] 
    for i in 1 : k
        push!(bases, x->sin(Ω*i*x[j]))
        push!(bases, x->cos(Ω*i*x[j])) 
    end
    return bases 
end
    
function sinusoidal_bases(k, Ω, a, b) 
    n = length(a)
    bases = [sinusoidal_bases_1d(i, k, Ω, a, b) for i in 1 : n] 
    terms = Function[]
    for ks in Iterators.product([0:2k for i in 1:n]...)
        powers = [div(k+1,2) for k in ks] 
        if sum(powers) ≤ k
            push!(terms, x->prod(b[j+1](x) for (j,b) in zip(ks,bases)))
        end 
    end
    return terms 
end


### Barke parameters 
@with_kw mutable struct p_barke 
    d :: Float64 = 100.e-9             # distance equilibrium cantilever to surface
    Q :: Float64 = 400.                # quality factor of the spring 
    k :: Float64 = 0.7                 # spring constant 
    ω_0 :: Float64 = 50.e3 * 2π        # eigenfrequency cantilever
    A_free :: Float64 = 250.e-9        # free amplitude
    force :: Float64 = k * A_free/Q
    Ω :: Float64 = 1.015               # frequency of the amplification
    V_0 :: Float64 = 4.2e-18           # Epsilon parameter bei Ingo, ACHTUNG er nutzt nur nm Skala --> 4.2 zu 4.2e-18
    σ :: Float64 = 2.8e-9              # Ingo (zz)
    δx :: Float64 = 0.5e-9             # softening parameter
    γ :: Float64 = 1.e-31              # 1.e-4 in Nm nm^2 s
    ϕ :: Float64 = 0.
end

@with_kw mutable struct p_DMT
    H :: Float64 = 2.e-19
    a_0 :: Float64 = 0.3e-9
    d :: Float64 = 8.5e-9
    E :: Float64 = 1.0905125408942204e9
    Q :: Float64 = 300.
    k :: Float64 = 40.
    Ω :: Float64 = 1.
    ω_0 :: Float64 = 2π * 300000.
    R :: Float64 = 10.e-9
    Γ :: Float64 = 1.240e-9/k            # = inputforce / k 
    ϕ :: Float64 = 0.
    i :: Int64 = 0
end

function affect!(integrator)
    # fourier series to gain the amplitude and phase from the last period
    integrator.p.i += 1
    base = sinusoidal_bases(2, integrator.p.Ω, integrator.t - t_period * Δt, integrator.t)
    t_b = collect(LinRange(integrator.t - t_period*Δt, integrator.t, floor(Int, t_period)))
    θ = lsq_regression(t_b, integrator.sol[1, end-floor(Int, t_period - 1.):end], base)
    A = sqrt(θ[2]^2 + θ[3]^2)
    φ = atan(θ[3], θ[2])
    # abspeichern der Ergebnisse 
    #push!(ampl_fun, A)
    ϕ = (integrator.t * integrator.p.Ω)%2π
    #push!(phi_fun, (ϕ - φ + π)%2π - π)
    popfirst!(std_buffer)
    push!(std_buffer, A)
    std_mean = std(std_buffer)/mean(std_buffer)    #push!(std_mean, std(std_buffer)/mean(std_buffer))
    #push!(std_mean_arr, std_mean)
    # push!(points_x, integrator.u[1])
    # push!(points_x_dot, integrator.u[2])
    # if integrator.t > t_min 
    #     # if norm([points_x[end], points_x_dot[end]] - [points_x[end-1], points_x_dot[end-1]]) < ϵ_tol
    #     #     terminate!(integrator)
    #     # end
    # end
    if log.(10, std_mean) < -5. && integrator.p.i > 30
        # if we converge for a given input force, we increase the the input force
        push!(phi_fun, (ϕ - φ + π)%2π - π)
        push!(ampl_fun, A)
        integrator.p.i = 0
        integrator.p.Γ += 8.163265306122441e-14
        #terminate!(integrator)
    end
    nothing
end

## DMT parameters
N = 50
f_high = 1.4e-9 
f_low = 1.240e-9 
f_fwd = collect(range(f_low, f_high, length = N))

#p = p_barke()
p = p_DMT()

### define staticarray problem
u0 = SA[0.0, 0.0]
tspan = (0., 100_000.)
Δt = 2π/1300
t_period = 1/p.Ω * 2π/Δt

ϵ_tol = 1e-16
ampl_fun = Float64[]
phi_fun = Float64[]

std_buffer = zeros(Float64, 30)
#std_mean_arr = Float64[]

periodic_cb = PeriodicCallback(affect!, t_period*Δt, initial_affect=false)
#stat_prob = ODEProblem(f_vLJ_static, u0, tspan, p)
stat_prob_DMT = ODEProblem(f_DMT, u0, tspan, p) 
integrator = init(stat_prob_DMT, callback=periodic_cb, alg=AutoTsit5(Rosenbrock23()), dt = Δt, adaptive=false)
solve!(integrator)



# ### some usual parameters for a simulation 
# p = [σ, δx, V_0, γ, 2π*f_0, Q, Ω, force/k, d, k, 0.]

# ### define staticarray problem
# u0 = SA[0.0, 0.0]
# tspan = (0., 8000.)
# Δt = 2π/1300
# t_period = 1/Ω * 2π/Δt
# std_buffer = zeros(Float64, 50)
# A_fin = Float64[]
# ϕ_fin = Float64[]
# ϕ = 0.
# sol_period = zeros(Float64, floor(Int, t_period + 1))


# function affect_only_period!(integrator)
#     base = sinusoidal_bases(2, Ω, integrator.t - t_period * Δt, integrator.t)
#     t_b = collect(LinRange(integrator.t - t_period*Δt, integrator.t, floor(Int, t_period)))
#     θ = lsq_regression(t_b, sol_period, base)
#     A = sqrt(θ[2]^2 + θ[3]^2)
#     ϕ = atan(θ[2]/θ[3])
#     popfirst!(std_buffer)
#     push!(std_buffer, A)
#     std_mean = std(std_buffer)/mean(std_buffer)
#     if log.(10, std_mean) < -4.
#         push!(A_fin, A)
#         push!(ϕ_fin, ϕ)
#         if length(A_fin) > 10
#             terminate!(integrator)
#         else
#             println(integrator.t)
#             old_ampl = integrator.p[8]
#             #reinit!(integrator, integrator.u ; t0 = integrator.t, tf = integrator.t + 3000.)
#             prob = remake(integrator.sol.prob; tspan=(integrator.t, integrator.t+3000), u0 = integrator.u)
#             integrator = init(prob, callback=periodic_cb, alg=AutoTsit5(Rosenbrock23()), dt = Δt, adaptive=false) 
#             integrator.p[8] = old_ampl + 0.05e-9
#         end
#     end
#     nothing
# end

# periodic_cb = PeriodicCallback(affect_only_period!, t_period*Δt, initial_affect=false)
# stat_prob = ODEProblem(f_vLJ_static, u0, tspan, p)
# integrator = init(stat_prob, callback=periodic_cb, alg=AutoTsit5(Rosenbrock23()), dt = Δt, adaptive=false)

# while integrator.t < 20000.
#     popfirst!(sol_period)
#     step!(integrator)
#     push!(sol_period, integrator.u[1])
# end





# A = ampl_fun[end]
# phi_test = phi_fun[end]
# fig = Figure(resolution=(800, 800))
# ax = Axis(fig[1, 1])
# t_test = collect(LinRange(integrator.t - t_period*Δt, integrator.t, floor(Int, t_period)))
# lines!(ax, t_test, -A.*cos.(Ω*t_test .- phi_test), linewidth = 3)
# hlines!(ax, ampl_fun[end])
# lines!(ax, t_test, integrator.sol[1, end-floor(Int, t_period - 1.):end], linestyle=:dash, linewidth=3)
# fig

# base = sinusoidal_bases(2, integrator.t - t_period * Δt, integrator.t)
# t_b = collect(integrator.t-t_period*Δt:Δt:integrator.t)
# θ = lsq_regression(t_b, integrator.sol[1, end-t_period:end], base)
