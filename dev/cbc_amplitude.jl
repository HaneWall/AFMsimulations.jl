using LinearAlgebra, CairoMakie

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

function sinusoidal_bases_1d(j, k, a, b) 
    T = b[j] - a[j]
    bases = Function[x->1/2] 
    for i in 1 : k
        push!(bases, x->sin(2π*i*x[j]/T))
        push!(bases, x->cos(2π*i*x[j]/T)) 
    end
    return bases 
end
    
function sinusoidal_bases(k, a, b) 
    n = length(a)
    bases = [sinusoidal_bases_1d(i, k, a, b) for i in 1 : n] 
    terms = Function[]
    for ks in Iterators.product([0:2k for i in 1:n]...)
        powers = [div(k+1,2) for k in ks] 
        if sum(powers) ≤ k
            push!(terms, x->prod(b[j+1](x) for (j,b) in zip(ks,bases)))
        end 
    end
    return terms 
end

function step_secant_reference(old::AbstractVector, μ::Float64)
    return 0 
end


function step_reference_trajectory(R::Float64, μ::Float64)
    R_new = R + μ
    return R_new
end

function fourier_decomposition(sig::AbstractVector; Ω_max::Int64 = 7)
    ampls_coeffs = lock_in_sig
    return coeff
end

function check_steady_state(fourier_coeffs_std::AbstractVector)
    return 0
end
