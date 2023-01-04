using DifferentialEquations, DiffEqCallbacks, StaticArrays, Parameters
using BenchmarkTools, LinearAlgebra, Statistics, CairoMakie, ProgressBars


""" 
returns effective Young's module of two interacting materials
"""
function eff_young_module(ν_t::Float64, E_t::Float64, ν_s::Float64, E_s::Float64)
    E = inv((1-ν_s^2)/E_s  + (1-ν_t^2)/E_t)
    return E
end

"""
This code should be able to reconstruct the results from Abeloos et. al. 
"""
@with_kw mutable struct p_dmt_ampl_sweeep 
    # Hamacker-Constant for Si-Si interaction
    H :: Float64 = 2.e-19
    # intermolecular distance
    a_0 :: Float64 = 0.3e-9
    # distance betweeen cantilever and surface/sample
    d :: Float64 = 8.5e-9 
    E_s :: Float64 = 1.e9
    ν_s :: Float64 = 0.3 
    E_t :: Float64 = 130.e9
    ν_t :: Float64 = 0.3
    # compute effective E with alrdy written function in utilities.jl 
    E :: Float64 = eff_young_module(ν_t, E_t, ν_s, E_s)
    # Quality-factor of cantilever
    Q :: Float64 = 300.
    # spring constant of cantilever
    k :: Float64 = 40.
    # resonance frequency
    f_0 :: Float64 = 300000.
    # Radius of Sphere
    R :: Float64 = 10.e-9
    # amplitude sweeps start and end 
    Γ_beg :: Float64 = 1. 
    Γ_end :: Float64 = 2.
    N :: Int64 = 80
    Γ_arr = collect(LinRange(Γ_beg, Γ_end, N))
    Γ :: Float64 = Γ_beg 
    # period counter in amplitude sweeps 
    i :: Int64 = 0
end


@with_kw mutable struct p_dmt_no_control 
    # Hamacker-Constant for Si-Si interaction
    H :: Float64 = 2.e-19
    # intermolecular distance
    a_0 :: Float64 = 0.3e-9
    # distance betweeen cantilever and surface/sample
    d :: Float64 = 8.5e-9 
    E_s :: Float64 = 1.e9
    ν_s :: Float64 = 0.3 
    E_t :: Float64 = 130.e9
    ν_t :: Float64 = 0.3
    # compute effective E with alrdy written function in utilities.jl 
    E :: Float64 = eff_young_module(ν_t, E_t, ν_s, E_s)
    # Quality-factor of cantilever
    Q :: Float64 = 300.
    # spring constant of cantilever
    k :: Float64 = 40.
    # resonance frequency
    f_0 :: Float64 = 300000.
    # Radius of Sphere
    R :: Float64 = 10.e-9
    
end

@with_kw mutable struct p_dmt_control
    # Hamacker-Constant for Si-Si interaction
    H :: Float64 = 2.e-19
    # intermolecular distance
    a_0 :: Float64 = 0.3e-9
    # distance betweeen cantilever and surface/sample
    d :: Float64 = 8.5e-9 
    E_s :: Float64 = 1.e9
    ν_s :: Float64 = 0.3 
    E_t :: Float64 = 130.e9
    ν_t :: Float64 = 0.3
    # compute effective E with alrdy written function in utilities.jl 
    E :: Float64 = eff_young_module(ν_t, E_t, ν_s, E_s)
    # Quality-factor of cantilever
    Q :: Float64 = 300.
    # spring constant of cantilever
    k :: Float64 = 40.
    # resonance frequency
    f_0 :: Float64 = 300000.
    # Radius of Sphere
    R :: Float64 = 10.e-9
    # GAIN-Parameters for the PD-controller 
    k_p :: Float64 = 0.
    k_d :: Float64 = 0.5
    # target-parameters for the control, X_1s_target = R_traget 
    X_1s_target :: Float64 = 1.
    X_1c_target :: Float64 = 0.
    # continuation parameter R_target_new = R_target + h = X_1c_target + h 
    h :: Float64 = 0.05
    # control input
    control :: Float64 = 0.
    N :: Int64 = 80
    # period counter in amplitude sweeps 
    i :: Int64 = 0
    # if we would like to record some values for debugging. This arrays are not
    # necessary if we found the right values. 
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
Functions that define our problem and the cbc-algorithm. 
"""

function f_dmt_no_control(u, p, t)
    @unpack k_1, k_2, k_3, Γ, Ω, ϕ = p
    dx = u[2]
    dy = - k_1 * u[1] - k_2 * u[2] - k_3 * u[1]^3 + Γ * sin(Ω*t + ϕ)
    return SA[dx, dy]
end


function f_dmt_no_control(x, p, t)
    @unpack Q, Ω, Γ, H, R, E, a_0, d, k_c, ϕ = p
    dx = x[2]
    dy = -1/Q * x[2] - x[1] + Γ*sin(Ω*t + ϕ)
    if x[1] < d - a_0
        dy += H*R/(6*k_c * (d - x[1])^2)
    else
        dy += H*R/(6 * a_0^2 * k_c) - 4/3*E*sqrt(R) * 1/(k_c) * (x[1] - (d - a_0))^(3/2) 
    end
    SA[dx, dy]
end 

function f_dmt(x, p, t)
    @unpack Q, control, H, R, E, a_0, d, k_c, ϕ = p
    dx = x[2]
    dy = -1/Q * x[2] - x[1] + control 
    if x[1] < d - a_0
        dy += H*R/(6*k_c * (d - x[1])^2)
    else
        dy += H*R/(6 * a_0^2 * k_c) - 4/3*E*sqrt(R) * 1/(k_c) * (x[1] - (d - a_0))^(3/2) 
    end
    SA[dx, dy]
end


function freq_sweep_stat(p::p_dmt_no_control)
    Δt = 0.005
    T_sim = 3000
    tspan = (0, T_sim)
    u_0 = SA[0.; 0.]
    ampl_container = zeros(Float64, p.N)
    phi_container = zeros(Float64, p.N)
    @inbounds for idx in ProgressBar(1:p.N)
        p.Ω = p.Ω_arr[idx]
        T_period = 2π/(Δt * p.Ω) 
        prob = ODEProblem(f_duffing_no_control, u_0, tspan, p)
        integrator = init(prob, AutoTsit5(Rosenbrock23()), dt=Δt, adaptive=false) 
        solve!(integrator)
        # amplitude detection algorithm (fft not necessary)
        T_period_int = ceil(Int, T_period)
        base = sinusoidal_bases(4, p.Ω, integrator.t - T_period * Δt, integrator.t)
        t_b = collect(LinRange(integrator.t - T_period*Δt, integrator.t, T_period_int))
        θ = lsq_regression(t_b, integrator.sol[1, end-(T_period_int-1):end], base)
        ampl_container[idx] = sqrt(θ[2]^2 + θ[3]^2)
        φ = atan(θ[3], θ[2]) 
        # phase difference input and output 
        phi_container[idx] = (p.ϕ - φ + π)%2π - π
        # phase conservation 
        p.ϕ = (p.ϕ + T_sim*p.Ω)%2π
        # new u_0 is last entry of previous simulation 
        u_0 = SA[integrator.sol[1, end]; integrator.sol[2, end]]
    end
    return ampl_container, phi_container
end

### amplitude sweeps, first no control, just up and down 
function  ampl_sweep(p::p_duffing_no_control_ampl_sweep)
    u0 =SA[0. ; 0.]
    Δt = 0.005
    tspan = (0., 25_000.)
    t_period = 1/p.Ω * 2π/Δt
    response_amplitude = Float64[]
    forcing_amplitude = Float64[]
   
    std_buffer = zeros(Float64, 30)
    function affect!(integrator)
        # this affect function takes place once every period (like a Poincare/stroboscopic map)
        # we count the periods i with the given parameters 
        integrator.p.i += 1
        # we get the amplitude and phase of the response from on the fly fourier coefficients
        A = sqrt(w[2]^2 + w[3]^2)
        # In order to detect the end of a transient behavior we create a 30
        # period buffer of the amplitude. If the standard deviation of this
        # buffer reaches a tolerance we converged. We can make this idea
        # dimensionless by division with the mean. 
        popfirst!(std_buffer)
        push!(std_buffer, A)
        std_mean = std(std_buffer)/mean(std_buffer)  
        if log.(10, std_mean) < -5.5 && integrator.p.i > 30
            # If we converge for a given input force, we increase the the input
            # force. 
            push!(response_amplitude, A)
            push!(forcing_amplitude, integrator.p.Γ)
            #restart period counter
            integrator.p.i = 0
            #continuation in excitation
            integrator.p.Γ += integrator.p.h
        end
        nothing
    end 

    periodic_cb = PeriodicCallback(affect!, t_period*Δt, initial_affect=false) # Poincare/Stroboscopic Map 
    prob = ODEProblem(f_duffing_no_control, u0, tspan, p) 
    integrator = init(prob, callback=periodic_cb, alg=AutoTsit5(Rosenbrock23()), dt = Δt, adaptive=false, save_everystep=false)

    k = 5                                           # approx signal up to 5th harmonic
    base = create_bases(k, integrator.p.Ω)          # basis for x online fourier
    
    # preallocate basis vectors and fourier coefficients
    w = zeros(2*k + 1)
    b = zeros(2*k + 1)
    steps = floor(Int64, tspan[2]/Δt)
    μ = Δt
    for st in ProgressBar(1:steps)
        update_basis!(b, base, integrator.t)
        update_weights!(w, b, integrator.u[1], μ)
        step!(integrator)
        if integrator.p.h > 0
            if integrator.p.Γ >= integrator.p.Γ_end
                terminate!(integrator)
                break
            end
        else
            if integrator.p.Γ <= integrator.p.Γ_end
                terminate!(integrator)
                break
            end 
        end
    end
    return response_amplitude, forcing_amplitude 
end

function cbc_sweep(p::p_abel_1)
    u0 = SA[0. ; 0.]
    Δt = 0.005
    tspan = (0., 25_000.)
    t_period = 1/p.Ω * 2π/Δt
    response_amplitude = Float64[]
    response_phase = Float64[]
    eff_forcing_amplitude = Float64[]

    t_buffer = zeros(Float64, ceil(Int64, t_period)) 
    Γ_buffer = zeros(Float64, ceil(Int64, t_period))
    std_buffer = zeros(Float64, 30)

    function affect!(integrator)
        # this affect function takes place every period (like a Poincare map)
        # we count the periods i with the given parameters 
        integrator.p.i += 1
        # we get the amplitude and phase of the response from on the fly fourier coefficients
        A = sqrt(w[2]^2 + w[3]^2)
        φ = atan(w[3], w[2])
        # We get the amplitude and phase of the effective forcing / control 
        ξ = lsq_regression(t_buffer, Γ_buffer, sin_base)
        Γ = sqrt(ξ[2]^2 + ξ[3]^2)
        ζ = (integrator.t * integrator.p.Ω)
        Γ_phi = (ζ - atan(ξ[3], ξ[2]) + π)%2π - π
        # In order to detect the end of a transient behavior we create a 30
        # period buffer of the amplitude. If the standard deviation of this
        # buffer reaches a tolerance we converged. We can make this idea
        # dimensionless by division with the mean. 
        popfirst!(std_buffer)
        push!(std_buffer, A)
        std_mean = std(std_buffer)/mean(std_buffer)  
        if log.(10, std_mean) < -5.5 && integrator.p.i > 30
            # If we converge for a given input force, we increase the the input
            # force. 
            push!(response_phase, φ)
            push!(response_amplitude, A)
            push!(eff_forcing_amplitude, Γ)
            #restart period counter
            integrator.p.i = 0
            #continuation in response
            integrator.p.X_1s_target += integrator.p.h
        end
        nothing
    end

    periodic_cb = PeriodicCallback(affect!, t_period*Δt, initial_affect=false) # Poincare/Stroboscopic Map 
    prob = ODEProblem(f_duffing, u0, tspan, p) 
    integrator = init(prob, callback=periodic_cb, alg=AutoTsit5(Rosenbrock23()), dt = Δt, adaptive=false, save_everystep=false)

    k = 5                                           # approx signal up to 5th harmonic
    base = create_bases(k, integrator.p.Ω)          # basis for x online fourier
    base_dot = create_dot_bases(k, integrator.p.Ω)  # basis for x_dot online fourier
    sin_base = sinusoidal_bases(k, integrator.p.Ω, 0, 2π) # fourier base - offline Least-Squares
    
    # preallocate basis vectors and fourier coefficients
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
        step!(integrator)
        if length(response_amplitude)>1 && eff_forcing_amplitude[end] >= 31.
            terminate!(integrator)
            break
        end
    end
    return response_amplitude, response_phase, eff_forcing_amplitude #, integrator.sol
end