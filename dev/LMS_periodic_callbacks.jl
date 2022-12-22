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
    @unpack Q, Ω, H, R, E, a_0, d, k, ϕ = p
    @unpack control = p
    dx = x[2]
    dy = -1/Q * x[2] - x[1] + 1/k * control
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

function p_controller(x_target::Float64, x::Float64, k_p::Float64)
    return k_p * (x_target - x)
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
        push!(bases, t-> -i * cos(Ω*i*t))
        push!(bases, t-> i * sin(Ω*i*t)) 
    end
    return bases
end

function update_basis!(b::AbstractArray{Float64}, bases::AbstractArray{Function}, t::Float64)
    b .= [f(t) for f in bases]
end


function update_weights!(w::AbstractArray{Float64}, b::AbstractArray{Float64}, data::Float64, μ::Float64)
    w .= w .+ μ .* inv(dot(b, b)) .* b .* ϵ_lms(data, w, b)
end

function update_controller(k_p::Float64, k_d::Float64, x_target::Float64, x::Float64, x_target_dot::Float64, x_dot::Float64)
    return k_p * (x_target - x) + k_d * (x_target_dot - x_dot)
end

"""
Following two functions are doing the same, but this might be better to
understand the underlying code.  
"""
function update_control_sin(k_p::Float64, X_target_1s::Float64, X_1s::Float64)
    return k_p * (X_target_1s - X_1s)
end

function update_control_cos(k_p::Float64, X_target_1c::Float64, X_1c::Float64)
    return k_p * (X_target_1c - X_1c)
end

"""
Γ is the new effective monoharmonic forcing term. --> needed for tracking. 
BTW: If there is no control, this is just the usual force. 
"""
function update_Γ!(Γ::Float64, force::Float64, control_sin::Float64, control_cos::Float64)
    Γ .= sqrt((force + control_sin)^2 + control_cos^2)
end

function calculate_x_nf(weights_nf::AbstractArray{Float64}, bases_nf::AbstractArray{Function}, t::Float64)
    b = [f(t) for f in bases_nf]
    return dot(weights_nf, b)
end

function update_x_target_nf!(x_target_nf::Float64, x_nf::Float64)
    x_target_nf .= x_nf
end

function continuation!(X_target_1s::Float64, h::Float64)
    X_target_1s .+= h
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

@with_kw mutable struct p_DMT_control 
    H :: Float64 = 2.e-19
    a_0 :: Float64 = 0.3e-9
    d :: Float64 = 8.5e-9
    E :: Float64 = 1.0905125408942204e9
    Q :: Float64 = 300.
    k :: Float64 = 40.
    Ω :: Float64 = 1.
    ω_0 :: Float64 = 2π * 300000.
    R :: Float64 = 10.e-9
    N :: Int64 = 80
    X_1s_target :: Float64 = 8.2e-9
    X_1c_target :: Float64 = 0
    h :: Float64 = 0.05e-9
    ϕ :: Float64 = 0.
    i :: Int64 = 0
    k_p :: Float64 = 60.
    k_d :: Float64 = 40.
    control :: Float64 = 0.
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
    for _ in ProgressBar(1:steps)
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


function amplitude_sweep_control(p::p_DMT_control)
    ### define staticarray problem
    u0 = SA[0.0, 0.0]
    tspan = (0., 60_000.)
    Δt = 2π/1300
    t_period = 1/p.Ω * 2π/Δt
    ampl_fun = Float64[]
    sizehint!(ampl_fun, 80)
    phi_fun = Float64[]
    sizehint!(phi_fun, 80) 
    Γ_fun = Float64[]
    sizehint!(Γ_fun, 80)
    ctrl = Float64[]
    std_buffer = zeros(Float64, 30)
    
    function affect!(integrator)
        integrator.p.i += 1
        A = sqrt(w[2]^2 + w[3]^2)
        φ = atan(w[3], w[2])
        ϕ = (integrator.t * integrator.p.Ω)%2π
        Γ = sqrt((p.k_p * (p.X_1s_target - w[2]*b[2]) - integrator.p.k_d*(p.X_1c_target - w[3]*b[3]))^2 + 
                 (p.k_p * (p.X_1c_target - w[3]*b[3]) + integrator.p.k_d*(p.X_1s_target - w[2]*b[2]))^2)
        popfirst!(std_buffer)
        push!(std_buffer, A)
        std_mean = std(std_buffer)/mean(std_buffer)  
        if log.(10, std_mean) < -5. && integrator.p.i > 30
            # if we converge for a given input force, we increase the the input force
            push!(phi_fun, (ϕ - φ + π)%2π - π)
            push!(ampl_fun, A)
            push!(Γ_fun, Γ)
            push!(ctrl, integrator.p.control)
            integrator.p.i = 0
            integrator.p.X_1s_target += integrator.p.h
        end
        nothing
    end
    periodic_cb = PeriodicCallback(affect!, t_period*Δt, initial_affect=false)
    stat_prob_DMT = ODEProblem(f_DMT_control, u0, tspan, p) 
    integrator = init(stat_prob_DMT, callback=periodic_cb, alg=AutoTsit5(Rosenbrock23()), dt = Δt, adaptive=false, save_everystep=false)

    k = 5
    base = create_bases(k, integrator.p.Ω)
    base_dot = create_dot_bases(k, integrator.p.Ω)
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
    for _ in ProgressBar(1:steps)
        step!(integrator)
        update_basis!(b, base, integrator.t)
        update_basis!(b_d, base_dot, integrator.t)
        update_weights!(w, b, integrator.u[1], Δt)
        update_nf_basis!(b_nf, b)
        update_nf_basis!(b_d_nf, b_d)
        update_nf_weights!(w_nf, w)
        x_nf = dot(w_nf, b_nf)
        x_target = x_nf + (integrator.p.X_1s_target*b[2] + integrator.p.X_1c_target*b[3])
        x_target_dot = integrator.p.Ω .* (integrator.p.X_1s_target*b_d[2] + integrator.p.X_1c_target*b_d[3]) + integrator.p.Ω * (dot(b_d_nf, w_nf))
        integrator.p.control = update_controller(integrator.p.k_p, integrator.p.k_d, x_target, integrator.u[1], x_target_dot, integrator.u[2])
        if length(ampl_fun) == integrator.p.N
            terminate!(integrator)
            break
        end
    end
    return ampl_fun, phi_fun, Γ_fun
end

function update_nf_basis!(b_nf::AbstractArray{Float64}, b::AbstractArray{Float64})
    b_nf .= [b[i] for i in eachindex(b) if i ∉ [2, 3]]
end 

function update_nf_weights!(w_nf::AbstractArray{Float64}, w::AbstractArray{Float64})
    w_nf .= [w[i] for i in eachindex(w) if i ∉ [2, 3]]
end 

# p_b_fwd = p_DMT()
# amplitudes_fwd, phis_fwd, forcing_fwd = amplitude_sweep_lms(p_b_fwd)

# p_b_bwd = p_DMT(F_beg=1.45e-9, F_end=1.15e-9)
# amplitudes_bwd, phis_bwd, forcing_bwd = amplitude_sweep_lms(p_b_bwd)

p_ctrl = p_DMT_control(k_p = 10., k_d = 30., h = 0.05e-9, X_1s_target=8.2e-9)
amplitudes_ctrl, phis_ctrl, eff_Γ = amplitude_sweep_control(p_ctrl)

CairoMakie.activate!(type="svg")
fig = Figure(resolution = (800, 800))
ax = Axis(fig[1, 1])
scatter!(ax, forcing_fwd, amplitudes_fwd)
scatter!(ax, forcing_bwd, amplitudes_bwd, markersize=5.5,color=:red)
scatter!(ax, eff_Γ, amplitudes_ctrl)
ax2 = Axis(fig[2, 1])
scatter!(ax2, forcing_fwd, phis_fwd)
scatter!(ax2, forcing_bwd, phis_bwd, markersize=5.5,color=:red)
scatter!(ax2, eff_Γ, phis_ctrl)
fig







