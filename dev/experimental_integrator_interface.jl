using DifferentialEquations, CairoMakie, StaticArrays, Statistics, LinearAlgebra


#=
First we should sweep without any continuation to gain information 
=#

k_max = 20 
ϕ̂_low = 0.3π
ϕ̂_high = 0.6π
ϕ̂_step = 0.1 
ϕ̂ = collect(ϕ̂_low:ϕ̂_step:ϕ̂_high)
K_p = 
K_d =  
j = 0
max_iter = 15

for ϕ_target in ϕ̂
    ϕ_exp = ϕ_exp
    ϕ̇_exp = ϕ̇_exp
    iter = 0
    while (norm(ϕ_exp - ϕ_target) > ϵ_ϕ || norm(ϕ̇_exp) > ϵ_ϕ̇) && iter < max_iter 
        # simulate transient time (10% of Q)
        integrator = init(ode_prob, )
        converged = false 
        while !converged
            step!(integrator)
        end
        # get current ϕ_exp, ϕ̇_exp, A_exp via lsq-regression?
        # increase ω
        ω = ω + T_trans * (K_p*(ϕ_exp - ϕ_target) + K_d*(ϕ̇_exp))
        iter += 1
    end
end


function ampl_sweep_init_scheme(F::Array{Float64}, Ω::Float64, exp::AFM_vLJ_experiment, Δt::Float64)
    N = length(F)
    ϕ = 0. 
    u_0 = SA[0.; 0.]
    ampl_container = zeros(Float64, N)
    Γ = F./(exp.tip.k) 
    @inbounds for idx in ProgressBar(1:N)
        p = [exp.σ, exp.δx, exp.V_0, exp.γ, 2π*exp.tip.f_0, exp.tip.Q, Ω, Γ[idx], exp.d, exp.tip.k, ϕ]
        prob = ODEProblem(f_vLJ!, u_0, tspan, p)
        cb = PeriodicCallback()
        integrator = init(prob, AutoTsit5(Rosenbrock23()), dt = Δt, adaptive=false)
        converged = false
        while !converged
            if converged
                terminate!(integrator)
                break 
            step!(integrator)
        end    
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