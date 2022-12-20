using LinearAlgebra, CairoMakie, DifferentialEquations, DiffEqCallbacks

"""
Control based continuation for the amplitude. 
Relevant paper:

Barton/Sieber:
https://arxiv.org/pdf/1209.3713

Springer, optimized:
https://link.springer.com/content/pdf/10.1007/s11071-021-06506-z.pdf?pdf=button

Barton:
https://arxiv.org/pdf/1506.04052.pdf
"""
Ω = 1. 
Δt = 0.005 
L = ceil(Int64, 2π * Ω / Δt) #length of one period
no_of_periods = 5
sig_container = zeros(Float64, no_of_periods * L)

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

function step_secant_reference(ξ, μ::Float64)
    """
    AbstractMatrix is of the form 
    ξ = [ [p_{n-1} x_{n-1}],  [p_{n}, x_{n}] ]
    """
    ξ_new = ξ[2] .+ μ.*(ξ[2] - ξ[1]) 
    popfirst!(ξ)
    push!(ξ, ξ_new)
    return ξ
end

function step_reference_trajectory(R::Float64, μ::Float64)
    R_new = R + μ
    return R_new
end

function fourier_decomposition(sig::AbstractVector; Ω_max::Int64 = 7)
    ampls_coeffs = lock_in_sig
    return coeff
end

"""
check every 10 periods, if we can increase the continuation. Callback on current
std_array --> if std_array[end] <= tol  
"""


function 

function check_steady_state(fourier_coeffs_std::AbstractVector)
    return 0
end



##PseudoCode for continuation in the amplitude: 
ampl_init = 1.                  # initial amplitude of the input signal 
N_H = 4                         # number of harmonic modes our basis spans, base = [1/2, sin(t), cos(t), .., sin(N_H*t), cos(N_H*t)] 
# x̂_ref = [X̂_0,X̂_1s ,X̂_1c , ... , X̂_NH_s, X̂_NH_c]    # reference trajectory
# x̂_nf_ref = [X̂_0, ... , X̂_NH_s, X̂_NH_c]  
# initial reference trajectory 
x̂_ref = zeros(Float64, 1+2*N_H)
x̂_ref[1] = ampl_init

u_0 = [0., 0.]                  # some initial condition 
integrator = ODEInit..          # setup an integrator that we can interate over
μ = 0.1                         # parameter for continuation for input amplitude 
tst_period = Ω*2π/Δt
sol_period = [zeros(Float64, tst_period), zeros(Float64, tst_period)]
std_buffer = zeros(Float64, 5)  # stores first harmonic for 5 period
ϵ_tol = 10e-10                  # Tolerance of standard deviation 
for tst in timesteps: 
    step!(integrator)
    popfirst!(sol_period)
    push!(sol_period, integrator.sol[end])
    # notice that we should use periodic saving Callbacks. (PeriodicCallback from DiffEqCallbacs)
    if tst % tst_period == 0
        base = sinusoidal_bases(freq_modes, (tst-tst_period)*Δt, tst*Δt)
        f_coeffs = lsq_regression(collect(tst-tst_period:1:tst).*Δt, sol_period, base)
        X_s, X_c = f_coeffs[2], f_coeffs[3]
        ϕ = atan(X_c/X_s)
        push!(fourier_coeffs, f_coeff)
        popfirst!(fourier_coeffs)
    elseif tst % (10*tst_period) == 0
        if check_steady_state(std_buffer)
            #- terminate integration 
            #- store current phase
            step_secant_reference(ξ, μ)
            break 
        else 
            continue 
        end
    end
end

mutable struct Reference_trajectory
    X̂ = zeros(Float64, Float64)

end


#computations that have to be done after each period
period_callback = PeriodicCallback(period_computations, timesteps_period; initialize=period_init)
period_compuations(integrator)



cbs = CallbackSet(period_callback, ten_period_callback)