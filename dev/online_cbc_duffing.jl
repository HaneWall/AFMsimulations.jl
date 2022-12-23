using DifferentialEquations, DiffEqCallbacks, StaticArrays, Parameters
using BenchmarkTools, LinearAlgebra, Statistics, CairoMakie, ProgressBars
 
"""
This code should be able to reconstruct some results from alrdy 
published cbc-experiments. 
"""

@with_kw mutable struct p_abel_1
    k_1 ::Float64 = 1.
    k_2 :: Float64 = 1. 
    k_3 :: Float64 = 1.
    k_p :: Float64 = 0.
    k_d :: Float64 = 0.5
    Ω :: Float64 = 4. 
    X_1s_target :: Float64 = 1.
    X_1c_target :: Float64 = 0. 
    h :: Float64 = 0.05
    control :: Float64 = 0.
    N :: Int64 = 80
    i :: Int64 = 0
    control_buff :: Vector{Float64} = Float64[]
    x_rec :: Vector{Float64} = Float64[]
    x_rec_dot :: Vector{Float64} = Float64[]
end

"""
Functions that enable on the fly fourier coefficients via adaptive notch filter.
"""

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

function ϵ_lms(data::Float64, weights::AbstractArray{Float64}, b::AbstractArray{Float64})
    return data - dot(weights, b)
end

function create_bases(k, Ω)
    bases = Function[t->1/2]
    for i in 1 : k
        push!(bases, t->sin(Ω*i*t))
        push!(bases, t->cos(Ω*i*t)) 
    end
    return bases
end

function create_dot_bases(k, Ω)
    bases = Function[t->0]
    for i in 1 : k
        push!(bases, t-> i * Ω * cos(Ω*i*t))
        push!(bases, t-> -i * Ω * sin(Ω*i*t)) 
    end
    return bases
end

function update_basis!(b::AbstractArray{Float64}, bases::AbstractArray{Function}, t::Float64)
    b .= [f(t) for f in bases]
end

function update_weights!(w::AbstractArray{Float64}, b::AbstractArray{Float64}, data::Float64, μ::Float64)
    w .= w .+ μ .* inv(dot(b, b)) .* b .* ϵ_lms(data, w, b)
end

function update_nf_basis!(b_nf::AbstractArray{Float64}, b::AbstractArray{Float64})
    b_nf .= [b[i] for i in eachindex(b) if i ∉ [2, 3]]
end 

function update_nf_weights!(w_nf::AbstractArray{Float64}, w::AbstractArray{Float64})
    w_nf .= [w[i] for i in eachindex(w) if i ∉ [2, 3]]
end 

"""
Functions that define ourt problem and the cbc-algorithm. 
"""

function f_duffing(u, p, t)
    @unpack k_1, k_2, k_3, control = p
    dx = u[2]
    dy = - k_1 * u[1] - k_2 * u[2] - k_3 * u[1]^3 + control 
    return SA[dx, dy]
end

function cbc_sweep(p::p_abel_1)
    u0 = SA[0. ; 0.]
    Δt = 0.005
    tspan = (0., 19_000.)
    t_period = 1/p.Ω * 2π/Δt
    response_amplitude = Float64[]
    response_phase = Float64[]
    eff_forcing_amplitude = Float64[]

    t_buffer = zeros(Float64, ceil(Int64, t_period)) 
    Γ_buffer = zeros(Float64, ceil(Int64, t_period))
    std_buffer = zeros(Float64, 30)
    
    function affect!(integrator)
        integrator.p.i += 1
        A = sqrt(w[2]^2 + w[3]^2)
        φ = atan(w[3], w[2])
        ϕ = (integrator.t * integrator.p.Ω)%2π
        #Γ = sqrt((integrator.p.k_p * (integrator.p.X_1s_target - w[2]*b[2]) - integrator.p.Ω * integrator.p.k_d*(integrator.p.X_1c_target - w[3]*b[3]))^2 + 
        #         (integrator.p.k_p * (integrator.p.X_1c_target - w[3]*b[3]) + integrator.p.Ω * integrator.p.k_d*(integrator.p.X_1s_target - w[2]*b[2]))^2)
        ξ = lsq_regression(t_buffer, Γ_buffer, sin_base)
        Γ = sqrt(ξ[2]^2 + ξ[3]^2)
        popfirst!(std_buffer)
        push!(std_buffer, A)
        std_mean = std(std_buffer)/mean(std_buffer)  
        if log.(10, std_mean) < -4.5 && integrator.p.i > 30
            # if we converge for a given input force, we increase the the input force
            push!(response_phase, (ϕ - φ + π)%2π - π)
            push!(response_amplitude, A)
            push!(eff_forcing_amplitude, Γ)
            integrator.p.i = 0
            integrator.p.X_1s_target += integrator.p.h
        end
        nothing
    end

    periodic_cb = PeriodicCallback(affect!, t_period*Δt, initial_affect=false)
    prob = ODEProblem(f_duffing, u0, tspan, p) 
    integrator = init(prob, callback=periodic_cb, alg=AutoTsit5(Rosenbrock23()), dt = Δt, adaptive=false, save_everystep=false)

    k = 5
    base = create_bases(k, integrator.p.Ω)
    base_dot = create_dot_bases(k, integrator.p.Ω)
    sin_base = sinusoidal_bases(k, integrator.p.Ω, 0, 2π)
    ξ = zeros(2*k + 1)
    w = zeros(2*k + 1)
    b = zeros(2*k + 1)
    b_d = zeros(2*k + 1)
    w_nf = zeros(2*k + 1 - 2)   # weights without fundamental
    b_nf = zeros(2*k + 1 - 2)   # basis without fundamental
    b_d_nf = zeros(2*k + 1 - 2)
    steps = floor(Int64, tspan[2]/Δt)
    x_nf = 0.
    x_target = 0.
    x_target_dot = 0.
    μ = Δt
    for st in ProgressBar(1:steps)
        update_basis!(b, base, integrator.t)
        update_basis!(b_d, base_dot, integrator.t)
        update_weights!(w, b, integrator.u[1], μ)

        update_nf_basis!(b_nf, b)
        update_nf_basis!(b_d_nf, b_d)
        update_nf_weights!(w_nf, w)
        x_target = p.X_1s_target*b[2] + dot(w_nf, b_nf)
        x_target_dot = p.X_1s_target * b_d[2] + dot(w_nf, b_d_nf)
        integrator.p.control =  (p.k_p * (x_target - integrator.u[1])) + p.k_d * (x_target_dot - integrator.u[2]) * tanh(st/(55*t_period)) 
        popfirst!(t_buffer)
        popfirst!(Γ_buffer)
        push!(t_buffer, integrator.t)
        push!(Γ_buffer, integrator.p.control)
        #push!(p.control_buff, integrator.p.control)
        step!(integrator)
        if length(response_amplitude)>1 && eff_forcing_amplitude[end] >= 31.
            terminate!(integrator)
            break
        end
    end
    return response_amplitude, response_phase, eff_forcing_amplitude, integrator.sol
end

omegas = [2. + i*0.25 for i in 1:1:8]
paras = [p_abel_1(X_1s_target=1., h=0.05, k_d=30.1, k_p=22.6, Ω=omegas[i]) for i in 1:8] 

#R_1, φ_R_1, Γ_1, sol_1 = cbc_sweep(paras[1])
R_2, φ_R_2, Γ_2 = cbc_sweep(paras[1])
R_3, φ_R_3, Γ_3 = cbc_sweep(paras[2])
R_4, φ_R_4, Γ_4 = cbc_sweep(paras[3])
R_5, φ_R_5, Γ_5 = cbc_sweep(paras[4])
R_6, φ_R_6, Γ_6 = cbc_sweep(paras[5])
R_7, φ_R_7, Γ_7 = cbc_sweep(paras[6])
R_8, φ_R_8, Γ_8 = cbc_sweep(paras[7])
R_9, φ_R_9, Γ_9 = cbc_sweep(paras[8])

CairoMakie.activate!(type="svg")
fig = Figure(resolution=(800, 800))
ax1 = Axis3(fig[1,1], azimuth=7.3*π/4, perspectiveness=0.)
#scatterlines!(ax1, Γ_1, R_1)
scatterlines!(ax1, 2.25 .* ones(length(R_2)), Γ_2, R_2)
scatterlines!(ax1, 2.5 .* ones(length(R_3)), Γ_3, R_3)
scatterlines!(ax1, 2.75 .* ones(length(R_4)), Γ_4, R_4)
scatterlines!(ax1, 3. .* ones(length(R_5)), Γ_5, R_5)
scatterlines!(ax1, 3.25 .* ones(length(R_6)), Γ_6, R_6)
scatterlines!(ax1, 3.5 .* ones(length(R_7)), Γ_7, R_7)
scatterlines!(ax1, 3.75 .* ones(length(R_8)), Γ_8, R_8)
scatterlines!(ax1, 4. .* ones(length(R_9)), Γ_9, R_9)
fig


# using GaussianProcesses
# using Plots
# ws = [2.5 .* ones(length(R_2)); 3. .* ones(length(R_3)); 3.5 .* ones(length(R_4)); 4. .* ones(length(R_5))]  
# Gs = [Γ_2; Γ_3; Γ_4; Γ_5]
# Rs = [R_2; R_3; R_4; R_5]
# Xs = transpose(hcat(Gs, ws))
# gp2 = GP(Xs,Rs, MeanZero(), SE(0., 0.))
# Plots.surface(gp2)

## testing dot base 

Ω = 7.3
x = 0.0:0.002:13π
y = 3. *sin.(Ω .* collect(x))  
y′ = 3. * Ω *cos.(Ω .* collect(x)) 

y_reconst = similar(x)
y_prime = similar(x)

k = 5
base = create_bases(k, Ω)
base_dot = create_dot_bases(k, Ω)
b = zeros(2*k + 1)
b_d = zeros(2*k + 1)
w = zeros(2*k + 1)

#
for (idx, t) in enumerate(x)
    update_basis!(b, base, t)
    update_basis!(b_d, base_dot, t)
    update_weights!(w, b, y[idx], 0.005)
    y_reconst[idx] = dot(w, b)
    y_prime[idx] = dot(w, b_d)
end

fig = Figure(resolution=(500, 500))
ax = Axis(fig[1, 1])
lines!(ax, x, y, color=:black)
lines!(ax, x, y′, color=:red)
lines!(ax, x, y_reconst, color=:black, linestyle=:dash)
lines!(ax, x, y_prime, color=:red, linestyle=:dash)
fig

function perfect_μ(Ω, Δt)
    return 2*sin(Ω*Δt) * (1 - sin(Ω * Δt))/(cos(Ω*Δt))^2
end