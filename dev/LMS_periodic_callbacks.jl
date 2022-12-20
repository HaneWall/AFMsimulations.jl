using DifferentialEquations, DiffEqCallbacks, StaticArrays, Parameters
using BenchmarkTools, LinearAlgebra, Statistics, CairoMakie, ProgressBars

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
    @unpack ξ_0, ξ_1s, ξ_1c, ξ_2s, ξ_2c, ξ_3s, ξ_3c = p  # target coeffs
    @unpack x_0, x_1s, x_1c, x_2s, x_2c, x_3s, x_3c = p  # current coeffs
    @unpack k_p, k_d = p                                 # gain factors
    dx = x[2]
    dy = -1/Q * x[2] - x[1] 
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
    F_beg :: Float64 = 1.15e-9
    F_end :: Float64 = 1.45e-9
    F :: Float64 = F_beg
    N :: Int64 = 80
    F_step = (F_end - F_beg)/N
    ϕ :: Float64 = 0.
    i :: Int64 = 0
end

# function affect!(integrator)
#     integrator.p.i += 1
#     A = sqrt(w[2]^2 + w[3]^2)
#     φ = atan(w[3], w[2])
#     ϕ = (integrator.t * integrator.p.Ω)%2π
#     popfirst!(std_buffer)
#     push!(std_buffer, A)
#     std_mean = std(std_buffer)/mean(std_buffer)    #push!(std_mean, std(std_buffer)/mean(std_buffer))
#     if log.(10, std_mean) < -5. && integrator.p.i > 30
#         # if we converge for a given input force, we increase the the input force
#         push!(phi_fun, (ϕ - φ + π)%2π - π)
#         push!(ampl_fun, A)
#         integrator.p.i = 0
#         integrator.p.Γ += integrator.p.Γ
#     end
#     nothing
# end

# ## DMT parameters
# N = 50
# f_high = 1.4e-9 
# f_low = 1.240e-9 
# f_fwd = collect(range(f_low, f_high, length = N))

# #p = p_barke()
# p = p_DMT(Γ = 1.240e-9/40.)

# ### define staticarray problem
# u0 = SA[0.0, 0.0]
# tspan = (0., 100_000.)
# Δt = 2π/1300
# t_period = 1/p.Ω * 2π/Δt

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