using DifferentialEquations
using ProgressBars
using NLsolve

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
    T_sim = 5000
    tspan = (0, T_sim)
    p = zeros(Float64, 10)
    u_0 = [0.; 0.]
    ampl_container = zeros(Float64, N)
    Γ = F/(exp.tip.k)
    @inbounds for idx in ProgressBar(1:N)
        p = [exp.tip.Q, Ω[idx], Γ, exp.sample.H, exp.tip.R, exp.eff_young, exp.a_0, exp.d, exp.tip.k, ϕ]
        prob = ODEProblem(f_DMT!, u_0, tspan, p)
        sol = solve(prob, AutoTsit5(Rosenbrock23()), dt=Δt, adaptive=false) 
        # amplitude detection algorithm (fft not necessary)
        T_period = Ω[idx] * 2π / Δt
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



"""
Frequency sweep of LJ-potential with fixpoints x^{⋆} of Poincare map Φ: x(t) -->
x(t + T) with T being the period of the forcing term. Should give even more
accurate results and might be more faster (due to less iterations).
"""
function freq_sweep_Φ(F::Float64, Ω::Array{Float64}, exp::AFM_vLJ_experiment, Δt::Float64)
    N = length(Ω)
    ϕ = 0. 
    u_0 = [0.; 0.]
    ampl_container = zeros(Float64, N)
    Γ = F/(exp.tip.k)
    @inbounds for idx in ProgressBar(1:N)
        p = [exp.σ, exp.δx, exp.V_0, exp.γ, 2π*exp.tip.f_0, exp.tip.Q, Ω[idx], Γ, exp.d, exp.tip.k, ϕ]
        T_period = Ω[idx] * 2π    
        tspan = (0. , 200. * T_period)
        T_sim = 200. * T_period
        prob_init = ODEProblem(f_vLJ!, u_0, tspan, p)
        sol_init = solve(prob_init, AutoTsit5(Rosenbrock23()), dt=Δt, adaptive=false) 
        #now we have some time integrated, however it might not been enough
        #for convergence. This is no problem tho, as long as we are in the basin
        #of attraction of the wanted fixpoint. We can now find the periodic
        #orbit, by finding the fixpoint of the Poincare map Φ_poincare. That is
        #we search for the solution x^{⋆} that fullfils of ϕ(x) - x = 0. When we
        #found x^{⋆} we can use it as the initial value and integrate over one
        #period to gain the amplitude.
        problem = ODEProblem(f_vLJ!, sol_init[end], (0., T_period), p)
        solution = solve(problem, AutoTsit5(Rosenbrock23()), dt=Δt, adaptive=false)
        ϕ = (ϕ + T_sim*Ω[idx])%2π
        p = [exp.σ, exp.δx, exp.V_0, exp.γ, 2π*exp.tip.f_0, exp.tip.Q, Ω[idx], Γ, exp.d, exp.tip.k, ϕ] 
        # Solve fixpoint problem with nonlinear solver
        solution2 = nlsolve(x -> Φ_poincare(problem, x, p, Δt) - x, solution[end])
        #ϕ = (ϕ + T_period*Ω[idx])%2π
        xstar = solution2.zero
        p = [exp.σ, exp.δx, exp.V_0, exp.γ, 2π*exp.tip.f_0, exp.tip.Q, Ω[idx], Γ, exp.d, exp.tip.k, ϕ]  

        final_solution = solve(remake(problem, u0=xstar, tspan=[0.0, T_period]), AutoTsit5((Rosenbrock23())), dt=Δt, adaptive=false)
        min = abs.(minimum(final_solution[1, :]))
        max = abs.(maximum(final_solution[1, :])) 
        ampl_container[idx] = 1/2 * (min + max)
        u_0 = [final_solution[1, end]; final_solution[2, end]]
        ϕ = (ϕ + T_period*Ω[idx])%2π
    end
    return ampl_container
end