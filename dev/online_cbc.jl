using DifferentialEquations, DiffEqCallbacks, StaticArrays, Parameters
using BenchmarkTools, LinearAlgebra, Statistics, CairoMakie, ProgressBars


"""
DMT model with StaticArrays.
"""
function f_DMT(x, p, t)
    @unpack Q, Ω, F, H, R, E, a_0, d, k, ϕ = p
    dx = x[2]
    dy = -1/Q * x[2] - x[1] + 1/k * F*sin(Ω*t + ϕ)
    if x[1] < d - a_0
        dy += H*R/(6*k * (d - x[1])^2)
    else
        dy += H*R/(6 * a_0^2 * k) - 4/3*E*sqrt(R) * 1/(k) * (x[1] - (d - a_0))^(3/2) 
    end
    SA[dx, dy]
end 

function f_DMT_control(x, p, t)
    @unpack Q, Ω, Γ, H, R, E, a_0, d, k, ϕ = p
    @unpack xᵗ = p  # target x
    @unpack k_p, k_d = p                                 # gain factors
    dx = x[2]
    dy = -1/Q * x[2] - x[1] + k_p * (ξ_nf - x[1]) 
    #dy += 1/k * k_p * ( (k*Γ - x_1s)*sin(Ω*t) + (0. - x_1c)*cos(Ω*t) -)
    if x[1] < d - a_0
        dy += H*R/(6*k * (d - x[1])^2)
    else
        dy += H*R/(6 * a_0^2 * k) - 4/3*E*sqrt(R) * 1/(k) * (x[1] - (d - a_0))^(3/2) 
    end
    SA[dx, dy]
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

function update_basis!(b::AbstractArray{Float64}, bases::AbstractArray{Function}, t::Float64)
    b .= [f(t) for f in bases]
end


function update_weights!(w::AbstractArray{Float64}, b::AbstractArray{Float64}, data::Float64, μ::Float64)
    w .= w .+ μ .* inv(dot(b, b)) .* b .* ϵ_lms(data, w, b)
end

# parameters Hölscher 
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
    F_beg :: Float64 = 1.15e-9
    F_end :: Float64 = 1.45e-9
    F :: Float64 = F_beg
    N :: Int64 = 80
    F_step = (F_end - F_beg)/N
    ϕ :: Float64 = 0.
    i :: Int64 = 0
end

function amplitude_sweep_lms(p::p_DMT)
    ### define staticarray problem
    u0 = SA[0.0, 0.0]
    tspan = (0., 120_000.)
    Δt = 2π/1300
    t_period = 1/p.Ω * 2π/Δt
    ampl_fun = Float64[]
    sizehint!(ampl_fun, 80)
    phi_fun = Float64[]
    sizehint!(phi_fun, 80) 
    forcing_fun = Float64[]
    sizehint!(forcing_fun, 80)
    std_buffer = zeros(Float64, 30)
    
    function affect!(integrator)
        integrator.p.i += 1
        A = sqrt(w[2]^2 + w[3]^2)
        φ = atan(w[3], w[2])
        ϕ = (integrator.t * integrator.p.Ω)%2π
        popfirst!(std_buffer)
        push!(std_buffer, A)
        std_mean = std(std_buffer)/mean(std_buffer)  
        if log.(10, std_mean) < -5. && integrator.p.i > 30
            # if we converge for a given input force, we increase the the input force
            push!(phi_fun, (ϕ - φ + π)%2π - π)
            push!(ampl_fun, A)
            push!(forcing_fun, integrator.p.F)
            integrator.p.i = 0
            integrator.p.F += integrator.p.F_step
        end
        nothing
    end
    periodic_cb = PeriodicCallback(affect!, t_period*Δt, initial_affect=false)
    stat_prob_DMT = ODEProblem(f_DMT, u0, tspan, p) 
    integrator = init(stat_prob_DMT, callback=periodic_cb, alg=AutoTsit5(Rosenbrock23()), dt = Δt, adaptive=false, save_everystep=false)

    k = 3
    base = create_bases(k, integrator.p.Ω)
    w = zeros(2*k + 1)
    b = zeros(2*k + 1)
    steps = floor(Int64, tspan[2]/Δt)
    for i in ProgressBar(1:steps)
        step!(integrator)
        update_basis!(b, base, integrator.t)
        update_weights!(w, b, integrator.u[1], Δt)
        if length(ampl_fun) == integrator.p.N
            terminate!(integrator)
            break
        end
    end
    return ampl_fun, phi_fun, forcing_fun
end
p_b_fwd = p_DMT()
amplitudes_fwd, phis_fwd, forcing_fwd = amplitude_sweep_lms(p_b_fwd)

p_b_bwd = p_DMT(F_beg=1.45e-9, F_end=1.15e-9)
amplitudes_bwd, phis_bwd, forcing_bwd = amplitude_sweep_lms(p_b_bwd)

CairoMakie.activate!(type="svg")
fig = Figure(resolution = (800, 800))
ax = Axis(fig[1, 1])
scatter!(ax, forcing_fwd, amplitudes_fwd)
scatter!(ax, forcing_bwd, amplitudes_bwd, markersize=5.5,color=:red)
ax2 = Axis(fig[2, 1])
scatter!(ax2, forcing_fwd, phis_fwd)
scatter!(ax2, forcing_bwd, phis_bwd, markersize=5.5,color=:red)
fig