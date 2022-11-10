using DifferentialEquations
using ProgressBars

"""
Forward or backward frequency sweep with constant excitation amplitude for LJ-model.
"""
function freq_sweep(F::Float64, Ω::Array{Float64}, exp::AFM_LJ_experiment, Δt::Float64)
    N = length(Ω)
    ϕ = 0.
    T_sim = 6000
    tspan = (0, T_sim)
    T_period = 2π/Δt
    p = zeros(Float64, 10)
    u_0 = [0.; 0.]
    ampl_container = zeros(Float64, N)
    Γ = F/(exp.tip.k)
    @inbounds for idx in ProgressBar(1:N)
        p = [exp.σ,  exp.tip.Q, Ω[idx], Γ, exp.sample.H, exp.tip.R, exp.eff_young, exp.a_0, exp.d, exp.tip.k, ϕ]
        prob = ODEProblem(f_LJ!, u_0, tspan, p)
        sol = solve(prob, AutoTsit5(Rosenbrock23()), dt=Δt, adaptive=false) 
        # amplitude detection algorithm (fft not necessary)
        timeslot_to_detect_ampl = ceil(Int, 2*T_period)
        min = abs.(minimum(sol[1, end-timeslot_to_detect_ampl:end]))
        max = abs.(maximum(sol[1, end-timeslot_to_detect_ampl:end]))
        ampl_container[idx] = 1/2 * (min + max)
        # phase conservation 
        ϕ = (ϕ + T_sim*Ω[idx])%2π
        # new u_0 is last entry of previous simulation 
        u_0 = [sol[1, end]; sol[2, end]]
    end
    return ampl_container
end


"""
Forward or backward frequency sweep with constant excitation amplitude for vLJ-model.
"""
function freq_sweep(F::Float64, Ω::Array{Float64}, exp::AFM_vLJ_experiment, Δt::Float64)
    N = length(Ω)
    ϕ = 0.
    T_sim = 4000
    tspan = (0, T_sim)
    T_period = 2π/Δt
    p = zeros(Float64, 10)
    u_0 = [0.; 0.]
    ampl_container = zeros(Float64, N)
    Γ = F/(exp.tip.k)
    @inbounds for idx in ProgressBar(1:N)
        p = [exp.σ, exp.δx, exp.V_0, exp.γ, 2π*exp.tip.f_0, exp.tip.Q, Ω[idx], Γ, exp.d, exp.tip.k, ϕ]
        prob = ODEProblem(f_vLJ!, u_0, tspan, p)
        sol = solve(prob, AutoTsit5(Rosenbrock23()), dt=Δt, adaptive=false) 
        # amplitude detection algorithm (fft not necessary)
        timeslot_to_detect_ampl = ceil(Int, 2*T_period)
        min = abs.(minimum(sol[1, end-timeslot_to_detect_ampl:end]))
        max = abs.(maximum(sol[1, end-timeslot_to_detect_ampl:end]))
        ampl_container[idx] = 1/2 * (min + max)
        # phase conservation 
        ϕ = (ϕ + T_sim*Ω[idx])%2π
        # new u_0 is last entry of previous simulation 
        u_0 = [sol[1, end]; sol[2, end]]
    end
    return ampl_container
end


"""
Forward or backward frequency sweep with constant excitation amplitude for DMT-model.
"""
function freq_sweep(F::Float64, Ω::Array{Float64}, exp::AFM_DMT_experiment, Δt::Float64)
    N = length(Ω)
    ϕ = 0.
    T_sim = 4000
    tspan = (0, T_sim)
    T_period = 2π/Δt
    p = zeros(Float64, 10)
    u_0 = [0.; 0.]
    ampl_container = zeros(Float64, N)
    Γ = F/(exp.tip.k)
    @inbounds for idx in ProgressBar(1:N)
        p = [exp.tip.Q, Ω[idx], Γ, exp.sample.H, exp.tip.R, exp.eff_young, exp.a_0, exp.d, exp.tip.k, ϕ]
        prob = ODEProblem(f_DMT!, u_0, tspan, p)
        sol = solve(prob, AutoTsit5(Rosenbrock23()), dt=Δt, adaptive=false) 
        # amplitude detection algorithm (fft not necessary)
        timeslot_to_detect_ampl = ceil(Int, 2*T_period)
        min = abs.(minimum(sol[1, end-timeslot_to_detect_ampl:end]))
        max = abs.(maximum(sol[1, end-timeslot_to_detect_ampl:end]))
        ampl_container[idx] = 1/2 * (min + max)
        # phase conservation 
        ϕ = (ϕ + T_sim*Ω[idx])%2π
        # new u_0 is last entry of previous simulation 
        u_0 = [sol[1, end]; sol[2, end]]
    end
    return ampl_container
end


"""
same for eDMT models
"""
function freq_sweep(F::Float64, Ω::Array{Float64}, exp::AFM_eDMT_experiment, Δt::Float64)
    N = length(Ω)
    ϕ = 0.
    T_sim = 6000
    tspan = (0, T_sim)
    T_period = 2π/Δt
    p = zeros(Float64, 10)
    u_0 = [0.; 0.]
    ampl_container = zeros(Float64, N)
    Γ = F/(exp.tip.k)
    @inbounds for idx in ProgressBar(1:N)
        p = [exp.γ_0, exp.x_γ, 2π*exp.tip.f_0, exp.tip.Q, Ω[idx], Γ, exp.sample.H, exp.tip.R, exp.eff_young, exp.a_0, exp.d, exp.tip.k, ϕ]
        prob = ODEProblem(f_eDMT!, u_0, tspan, p)
        sol = solve(prob, AutoTsit5(Rosenbrock23()), dt=Δt, adaptive=false) 
        # amplitude detection algorithm (fft not necessary)
        timeslot_to_detect_ampl = ceil(Int, 2*T_period)
        min = abs.(minimum(sol[1, end-timeslot_to_detect_ampl:end]))
        max = abs.(maximum(sol[1, end-timeslot_to_detect_ampl:end]))
        ampl_container[idx] = 1/2 * (min + max)
        # phase conservation 
        ϕ = (ϕ + T_sim*Ω[idx])%2π
        # new u_0 is last entry of previous simulation 
        u_0 = [sol[1, end]; sol[2, end]]
    end
    return ampl_container
end


"""
same for vDMT models
"""
function freq_sweep(F::Float64, Ω::Array{Float64}, exp::AFM_vDMT_experiment, Δt::Float64)
    N = length(Ω)
    ϕ = 0.
    T_sim = 6000
    tspan = (0, T_sim)
    T_period = 2π/Δt
    p = zeros(Float64, 10)
    u_0 = [0.; 0.]
    ampl_container = zeros(Float64, N)
    Γ = F/(exp.tip.k)
    @inbounds for idx in ProgressBar(1:N)
        p = [exp.η, 2π*exp.tip.f_0, exp.tip.Q, Ω[idx], Γ, exp.sample.H, exp.tip.R, exp.eff_young, exp.a_0, exp.d, exp.tip.k, ϕ]
        prob = ODEProblem(f_eDMT!, u_0, tspan, p)
        sol = solve(prob, AutoTsit5(Rosenbrock23()), dt=Δt, adaptive=false) 
        # amplitude detection algorithm (fft not necessary)
        timeslot_to_detect_ampl = ceil(Int, 2*T_period)
        min = abs.(minimum(sol[1, end-timeslot_to_detect_ampl:end]))
        max = abs.(maximum(sol[1, end-timeslot_to_detect_ampl:end]))
        ampl_container[idx] = 1/2 * (min + max)
        # phase conservation 
        ϕ = (ϕ + T_sim*Ω[idx])%2π
        # new u_0 is last entry of previous simulation 
        u_0 = [sol[1, end]; sol[2, end]]
    end
    return ampl_container
end