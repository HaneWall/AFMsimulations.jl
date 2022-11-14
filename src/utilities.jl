using Statistics, DifferentialEquations, CairoMakie
CairoMakie.activate!(type="svg")

""" 
returns effective Young's module of two interacting materials
"""
function eff_young_module(s::Sample, t::Cantilever)
    E = ( (1-s.ν^2)/s.E  + (1-t.ν^2)/t.E )^(-1)
    return E
end


"""
returns force-distance relation for DMT experiments
"""
function force_distance(x::Float64, exp::AFM_DMT_experiment)
    R = exp.tip.R
    H = exp.sample.H
    E = exp.eff_young
    d = exp.d
    a_0 = exp.a_0
    if x < d - a_0
        return -H*R/(6*(d-x)^2)
    else 
        return -H*R/(6 * a_0^2) + 4/3*E*sqrt(R)* (x - (d - a_0))^(3/2) 
    end
end


"""
returns force distance relation for LJ experiments
"""
function force_distance(x::Float64, exp::AFM_vLJ_experiment)
    σ = exp.σ
    V_0 = exp.V_0 
    δx = exp.δx 
    h_x = x^2 + δx^2
    f = 12*V_0/sqrt(h_x) *( (σ^2/ h_x)^6 - (σ^2/ h_x)^3)
    return f
end

"""
Poincare map Φ: x(t) --> x(t + T), with T being the period of the forcing term.
For us T = 2π * Ω. Input is a ODE-Problem, initial Vector x, problem specific
parameters and Δt.
"""
function Φ_poincare(prob::ODEProblem, x::Array{Float64}, p::Array{Float64}, Δt::Float64)
    prob_new = remake(prob, u0=x, p=p)
    sol_new = solve(prob_new, AutoTsit5(Rosenbrock23()), adaptive=false, dt=Δt) 
    return sol_new[end]
end




"""
Steady state checker, stop simulation if std of amplitude for n periods stays
below some tolerance: tol  
"""
function steady_state_check(ampl::Array{Float64}, tol::Float64)
    σ = std(ampl)
    if σ < tol
        return true
    else 
        return false
    end
end


"""
Implements digital lock in amplifier. Provide a reference signal (a local
osciallator with ω=ω_ref) and gain information of the phase and amplitude of the
input signal. 
Input:: sig:Array of interest, ω_ref:reference frequency 
"""
function lock_in_amplifier(sig::Array{Float64}, ω_ref::Float64)
end 


"""
Creates phase space representation for the damped Lennard-Jones model.
"""
function phase_space(exp::AFM_vLJ_experiment)
    x_min = -exp.d
    x_max = exp.d
    ∂x_min = -10. *exp.d
    ∂x_max = 10. *exp.d 
    h_x(x) = (x - exp.d)^2 + (exp.δx)^2
    help_function(x, y) = Point2f(y, -1/exp.tip.Q * y - x - 12*exp.V_0/(exp.tip.k * sqrt(h_x(x))) * ((exp.σ^2 / h_x(x))^6 - (exp.σ^2 / h_x(x))^3)  - y * exp.tip.f_0*2π * exp.γ/(exp.d - x)^3)
    fig = Figure(resolution=(500, 500), title="Phase Space for damped Lennard-Jones model")
    ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"∂x")
    streamplot!(ax, help_function, x_min .. x_max, ∂x_min .. ∂x_max, colormap = [:black, :black], gridsize = (10, 10), arrow_size = 10)
    fig
end