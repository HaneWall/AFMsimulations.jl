using DifferentialEquations, DiffEqCallbacks, StaticArrays, Parameters
using BenchmarkTools, LinearAlgebra, Statistics, CairoMakie, ProgressBars
 

@with_kw mutable struct p_no_control_ampl_sweep
    k_1 :: Float64 = 20. / 20.
    k_2 :: Float64 = 0.5 / 20. 
    k_3 :: Float64 = 1. / 20.
    N :: Int64 = 80
    i :: Int64 = 0
    ϕ :: Float64 = 0. 
    Γ_beg :: Float64 = 0.5
    Γ_end :: Float64 = 23. 
    h :: Float64 = (Γ_end - Γ_beg)/N
    Γ :: Float64 = Γ_beg
    Ω :: Float64 = 0.87
end


@with_kw mutable struct p_no_control
    k_1 :: Float64 = 20. / 20. 
    k_2 :: Float64 = 0.5 / 20.
    k_3 :: Float64 = 1. / 20.
    Ω_beg :: Float64 = 0.7 
    Ω_end :: Float64 = 1.1
    N :: Int64 = 80
    i :: Int64 = 0
    ϕ :: Float64 = 0. 
    Γ :: Float64 = 10.
    Ω_arr :: Vector{Float64} = collect(range(Ω_beg, Ω_end, length=N))
    Ω :: Float64 = Ω_beg
end

@with_kw mutable struct p_control
    k_1 ::Float64 = 20. / 20.
    k_2 :: Float64 = 0.5 / 20. 
    k_3 :: Float64 = 1. / 20.
    k_p :: Float64 = 0.
    k_d :: Float64 = 3.
    Ω :: Float64 = 0.87 
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

function f_no_control(u, p, t)
    @unpack k_1, k_2, k_3, Γ, Ω, ϕ = p
    dx = u[2]
    dy = - k_1 * u[1] - k_2 * u[2] - k_3 * (10 * u[1]^2 + 2 * u[1]^3) + k_3 * Γ * sin(Ω*t + ϕ)
    return SA[dx, dy]
end


function f_control(u, p, t)
    @unpack k_1, k_2, k_3, control = p
    dx = u[2]
    dy = - k_1 * u[1] - k_2 * u[2] - k_3 * (10 * u[1]^2 + 2 * u[1]^3) + k_3 * control 
    return SA[dx, dy]
end


function freq_sweep_stat(p::p_no_control)
    Δt = 0.005
    T_sim = 3000
    tspan = (0, T_sim)
    u_0 = SA[0.; 0.]
    ampl_container = zeros(Float64, p.N)
    phi_container = zeros(Float64, p.N)
    @inbounds for idx in ProgressBar(1:p.N)
        p.Ω = p.Ω_arr[idx]
        T_period = 2π/(Δt * p.Ω) 
        prob = ODEProblem(f_no_control, u_0, tspan, p)
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
function  ampl_sweep(p::p_no_control_ampl_sweep)
    u0 =SA[0. ; 0.]
    Δt = 0.005
    tspan = (0., 55_000.)
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
        if log.(10, std_mean) < -4.5 && integrator.p.i > 30
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
    prob = ODEProblem(f_no_control, u0, tspan, p) 
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

function cbc_sweep(p::p_control)
    u0 = SA[0. ; 0.]
    Δt = 0.005
    tspan = (0., 50_000.)
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
        if log.(10, std_mean) < -4.5 && integrator.p.i > 30
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
    prob = ODEProblem(f_control, u0, tspan, p) 
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

## some amplitude sweeps without control 

p_ampl_087_fwd = p_no_control_ampl_sweep(Γ_beg = 0.2 , Γ_end = 7., N=100, Ω = 0.88 )
p_ampl_087_bwd = p_no_control_ampl_sweep(Γ_beg = 7., Γ_end = 0.2, N=100, Ω = 0.88 )

a_fwd, Γ_fwd = ampl_sweep(p_ampl_087_fwd)
a_bwd, Γ_bwd = ampl_sweep(p_ampl_087_bwd) 

# compare with cbc
p_cbc = p_control(X_1s_target = 0.1, Ω = 0.88, k_p=3., k_d=20.)

a_cbc, phi_cbc, Γ_cbc = cbc_sweep(p_cbc)

CairoMakie.activate!(type="svg")
fig = Figure(resolution = (800, 800), fontsize=28, fonts = (; regular = "CMU Serif"))
ax = Axis(fig[1,1], xlabel=L"\Gamma", ylabel=L"R")
scatter!(ax, Γ_fwd, a_fwd, marker=:rect, markersize=24, strokecolor=:black, strokewidth=2, color=(:black, 0.), label="Forward")
scatter!(ax, Γ_bwd, a_bwd, marker=:circle, markersize=24, strokecolor=:red, strokewidth=2, color=(:black, 0.), label="Backward")
scatterlines!(ax, Γ_cbc, a_cbc, label="CBC")
axislegend(L"\omega = 0.88", position=:rc)
xlims!(ax, (0.2, 4.))
ylims!(ax, (-0.1, 3.))
fig


## some frequency sweeps to compare theory to real sweeps
p_no_fwd = p_no_control(Ω_beg = 0.7, Ω_end = 1.1, Γ = 1.6)
p_no_bwd = p_no_control(Ω_end = 0.7 , Ω_beg = 1.1, Γ = 1.6)

r_fwd, phi_fwd = freq_sweep_stat(p_no_fwd)
r_bwd, phi_bwd = freq_sweep_stat(p_no_bwd)

CairoMakie.activate!(type="svg")
fig = Figure(resolution = (800, 800), fontsize=28, fonts = (; regular = "CMU Serif"))
ax = Axis(fig[1,1], xlabel=L"\omega", ylabel=L"R")
scatter!(ax, range(0.7, 1.1, 80), r_fwd, marker=:rect, markersize=24, strokecolor=:black, strokewidth=2, color=(:black, 0.), label="Forward")
scatter!(ax, range(1.1, 0.7, 80), r_bwd, marker=:circle, markersize=24, strokecolor=:red, strokewidth=2, color=(:black, 0.), label="Backward")
axislegend(L"\Gamma= 1.6", position=:rc)
fig

# p_no_fwd_9 = p_duffing_no_control(Γ = 9.)
# p_no_bwd_9 = p_duffing_no_control(Ω_beg = 4., Ω_end=1.5, Γ = 9.)

# p_no_fwd_45 = p_duffing_no_control(Γ = 4.5)
# p_no_bwd_45 = p_duffing_no_control(Ω_beg = 4., Ω_end=1.5, Γ = 4.5)



# r_fwd_9, phi_fwd_9 = freq_sweep_stat(p_no_fwd_9)
# r_bwd_9, phi_bwd_9 = freq_sweep_stat(p_no_bwd_9)

# r_fwd_45, phi_fwd_45 = freq_sweep_stat(p_no_fwd_45)
# r_bwd_45, phi_bwd_45 = freq_sweep_stat(p_no_bwd_45)


## amplitude sweeps with control 
omegas = collect(range(0.8, 1., length=20))
paras = [p_control(X_1s_target=0.1, h=0.05, k_d=8.1, k_p=0., Ω=omegas[i]) for i in 1:20] 

# TODO: change names, fked up with the numbers
R_2, φ_R_2, Γ_2 = cbc_sweep(paras[1])
R_3, φ_R_3, Γ_3 = cbc_sweep(paras[2])
R_4, φ_R_4, Γ_4 = cbc_sweep(paras[3])
R_5, φ_R_5, Γ_5 = cbc_sweep(paras[4])
R_6, φ_R_6, Γ_6 = cbc_sweep(paras[5])
R_7, φ_R_7, Γ_7 = cbc_sweep(paras[6])
R_8, φ_R_8, Γ_8 = cbc_sweep(paras[7])
R_9, φ_R_9, Γ_9 = cbc_sweep(paras[8])
R_10, φ_R_10, Γ_10 = cbc_sweep(paras[9])
R_11, φ_R_11, Γ_11 = cbc_sweep(paras[10])
R_12, φ_R_12, Γ_12 = cbc_sweep(paras[11])
R_13, φ_R_13, Γ_13 = cbc_sweep(paras[12])
R_14, φ_R_14, Γ_14 = cbc_sweep(paras[13])
R_15, φ_R_15, Γ_15 = cbc_sweep(paras[14])
R_16, φ_R_16, Γ_16 = cbc_sweep(paras[15])
R_17, φ_R_17, Γ_17 = cbc_sweep(paras[16])
R_18, φ_R_18, Γ_18 = cbc_sweep(paras[17])
R_19, φ_R_19, Γ_19 = cbc_sweep(paras[18])
R_20, φ_R_20, Γ_20 = cbc_sweep(paras[19])
R_21, φ_R_21, Γ_21 = cbc_sweep(paras[20])


using GaussianProcesses
ws = [omegas[1] .* ones(length(R_2)); omegas[2] * ones(length(R_3)); omegas[3] .* ones(length(R_4)); omegas[4] .* ones(length(R_5)); omegas[5] .* ones(length(R_6)); omegas[6] .* ones(length(R_7)); omegas[7] .* ones(length(R_8)); omegas[8] .* ones(length(R_9)); omegas[9] .* ones(length(R_10)); omegas[10] .* ones(length(R_11)); omegas[11] .* ones(length(R_12)); omegas[12] .* ones(length(R_13)); omegas[13] .* ones(length(R_14)); omegas[14] .* ones(length(R_15));omegas[15] .* ones(length(R_16)); omegas[16] .* ones(length(R_17));omegas[17] .* ones(length(R_18)); omegas[18] .* ones(length(R_19)); omegas[19] .* ones(length(R_20)); omegas[20] .* ones(length(R_21)) ]  
Gs = [Γ_2; Γ_3; Γ_4; Γ_5; Γ_6; Γ_7; Γ_8; Γ_9; Γ_10; Γ_11; Γ_12; Γ_13; Γ_14; Γ_15; Γ_16; Γ_17; Γ_18; Γ_19; Γ_20; Γ_21]
Rs = [R_2; R_3; R_4; R_5; R_6; R_7; R_8; R_9; R_10; R_11; R_12; R_13; R_14; R_15; R_16; R_17; R_18; R_19; R_20; R_21]
φs = [φ_R_2; φ_R_3; φ_R_4; φ_R_5; φ_R_6; φ_R_7; φ_R_8; φ_R_9; φ_R_10; φ_R_11; φ_R_12; φ_R_13; φ_R_14; φ_R_15; φ_R_16; φ_R_17; φ_R_18; φ_R_19; φ_R_20; φ_R_21] 

Xs = transpose(hcat(Rs, ws))
Xs_phi = transpose(hcat(φs, ws))
gp2 = GP(Xs,Gs, MeanZero(), SE(-2.8, 0.))
gpphi = GP(Xs_phi, Gs, MeanZero(), SE(0., 0.)) 

xmin, xmax =  (minimum(gp2.x[1,:]), maximum(gp2.x[1,:]))
ymin, ymax = (minimum(gp2.x[2,:]), maximum(gp2.x[2,:]))
x = range(xmin, stop=xmax, length=250)
y = range(ymin, stop=ymax, length=250)
xgrid = repeat(x', 250, 1)
ygrid = repeat(y, 1, 250)
μ, Σ = predict_f(gp2,[vec(xgrid)'; vec(ygrid)'])
zgrid = reshape(μ, 250, 250)


using GLMakie
GLMakie.activate!()
z_index = Observable(1.)
lvls = @lift(range(0., $z_index, step=$z_index))
fig = Figure(resolution=(1400, 700), fontsize=24, fonts = (; regular = "CMU Serif"))
#ax = Axis3(fig[1, 1], azimuth = 0., elevation=-π/2) 
ax = Axis3(fig[1, 1], xlabel=L"R", ylabel=L"\omega", zlabel=L"\Gamma")
CairoMakie.surface!(ax, xgrid, ygrid, zgrid, colormap=:Spectral_11, shading=false)
scatterlines!(ax, R_2, omegas[1] .* ones(length(R_2)), Γ_2, color=:black)
scatterlines!(ax, R_3, omegas[2] .* ones(length(R_3)), Γ_3, color=:black)
scatterlines!(ax, R_4, omegas[3] .* ones(length(R_4)), Γ_4, color=:black)
scatterlines!(ax, R_5, omegas[4] .* ones(length(R_5)), Γ_5, color=:black)
scatterlines!(ax, R_6, omegas[5] .* ones(length(R_6)), Γ_6, color=:black)
scatterlines!(ax, R_7, omegas[6] .* ones(length(R_7)), Γ_7, color=:black)
scatterlines!(ax, R_8, omegas[7] .* ones(length(R_8)), Γ_8, color=:black)
scatterlines!(ax, R_9, omegas[8] .* ones(length(R_9)), Γ_9, color=:black)
scatterlines!(ax, R_10, omegas[9] .* ones(length(R_10)), Γ_10, color=:black)
scatterlines!(ax, R_11, omegas[10] .* ones(length(R_11)), Γ_11, color=:black)
scatterlines!(ax, R_12, omegas[11] .* ones(length(R_12)), Γ_12, color=:black)
scatterlines!(ax, R_13, omegas[12] .* ones(length(R_13)), Γ_13, color=:black)
scatterlines!(ax, R_14, omegas[13] .* ones(length(R_14)), Γ_14, color=:black)
scatterlines!(ax, R_15, omegas[14] .* ones(length(R_15)), Γ_15, color=:black)
scatterlines!(ax, R_16, omegas[15] .* ones(length(R_16)), Γ_16, color=:black)
scatterlines!(ax, R_17, omegas[16] .* ones(length(R_17)), Γ_17, color=:black)
scatterlines!(ax, R_18, omegas[17] .* ones(length(R_18)), Γ_18, color=:black)
scatterlines!(ax, R_19, omegas[18] .* ones(length(R_19)), Γ_19, color=:black)
scatterlines!(ax, R_20, omegas[19] .* ones(length(R_20)), Γ_20, color=:black)
scatterlines!(ax, R_21, omegas[20] .* ones(length(R_21)), Γ_21, color=:black)
CairoMakie.contour3d!(ax, x, y, zgrid'; levels=lvls, color=:white, linewidth=5.)
ax2 = Axis(fig[1, 2], xlabel=L"\omega", ylabel=L"R")
GLMakie.contour!(ax2, y, x, zgrid, levels=lvls, color=:black, linewidth=3.)
sl = Slider(fig[1, 3], horizontal = false, range = 0.1:0.1:2.)
ylims!(ax2, (0.2, 3.))
connect!(z_index, sl.value)
fig


CairoMakie.activate!(type="svg")
fig_heat = Figure(resolution=(800, 700), fontsize=28, fonts = (; regular = "CMU Serif"))
ax2 = Axis(fig_heat[1, 1], xlabel=L"\omega", ylabel=L"R")
heater = CairoMakie.heatmap!(ax2, ygrid, xgrid, zgrid; levels=1.:0.2:2., colormap=Reverse(:Spectral_11), colorrange=(1., 2.), highclip=:gray) 
CairoMakie.contour!(ax2, y, x, zgrid; levels=1.:0.2:2., color=:black)
scatter!(ax2, p_no_bwd.Ω_arr, r_bwd, marker=:ltriangle, color=:red, label="Backward")
scatter!(ax2, p_no_fwd.Ω_arr, r_fwd, marker=:rtriangle, color=:black, label="Forward")
CairoMakie.xlims!(ax2, [0.8, 1.])
CairoMakie.ylims!(ax2, [0.2, 3.])
axislegend(position=:lt)
Colorbar(fig_heat[1, 2][1, 1], heater, ticks=1.:0.2:2., vertical=true)
Label(fig_heat[1, 2][1, 2], L"\Gamma", tellheight=false, rotation=π/2)
fig_heat

GLMakie.activate!()
fig_trans = Figure(resolution = (800, 800)) 
ax_trans = Axis3(fig_trans[1, 1])
surface!(ax_trans, xgrid, zgrid, ygrid, colormap=[:gray, :gray], transparency=true, shading=true)
ylims!(ax_trans, (0, 10.))
xlims!(ax_trans, (0., 3.5))
fig_trans



## experiments with otf fourier
# Ω = 7.3
# x = 0.0:0.002:13π
# y = 3. *sin.(Ω .* collect(x))  
# y′ = 3. * Ω *cos.(Ω .* collect(x)) 

# y_reconst = similar(x)
# y_prime = similar(x)

# k = 5
# base = create_bases(k, Ω)
# base_dot = create_dot_bases(k, Ω)
# b = zeros(2*k + 1)
# b_d = zeros(2*k + 1)
# w = zeros(2*k + 1)

# #
# for (idx, t) in enumerate(x)
#     update_basis!(b, base, t)
#     update_basis!(b_d, base_dot, t)
#     update_weights!(w, b, y[idx], 0.005)
#     y_reconst[idx] = dot(w, b)
#     y_prime[idx] = dot(w, b_d)
# end

# fig = Figure(resolution=(500, 500))
# ax = Axis(fig[1, 1])
# lines!(ax, x, y, color=:black)
# lines!(ax, x, y′, color=:red)
# lines!(ax, x, y_reconst, color=:black, linestyle=:dash)
# lines!(ax, x, y_prime, color=:red, linestyle=:dash)
# fig

# function perfect_μ(Ω, Δt)
#     return 2*sin(Ω*Δt) * (1 - sin(Ω * Δt))/(cos(Ω*Δt))^2
# end