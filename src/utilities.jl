using Statistics, DifferentialEquations, CairoMakie, LinearAlgebra, DSP
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
Poincare map Φ: x(t) --> x(t + T), with T being the period of the forcing term.
For us T = 2π * Ω. Input is a ODE-Problem, initial Vector x, problem specific
parameters and Δt. We found a periodic orbit iff |Φ(x) - x|_2 <= ϵ with ϵ being
the absolute tolerance. This function can tzhan be used to determine a callback
function in an ODEProblem. This way we Should be able to get way more efficient
algorithms. 
"""
function Φ_poincare_fixpoint(prob::ODEProblem, x::Array{Float64}, p::Array{Float64}, Δt::Float64, abstol::Float64)
    prob_new = remake(prob, u0=x, p=p)
    sol_new = solve(prob_new, AutoTsit5(Rosenbrock23()), adaptive=false, dt=Δt)  
    if norm(sol_new[end] .- x).< abstol
        return true
    else 
        return false
    end
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


"""
Savitzky-Golay Filter.
"""
function savitzky_golay_filter(y::AbstractVector, window_size::Integer, polynomial_order::Integer; deriv_order::Integer = 0, boundary_mode = :interpolation)

	# raise errors if input is bad
   	@assert isodd(window_size) "Window size must be an odd integer, i.e. fitting 2m + 1 points around the current value."
	@assert polynomial_order < window_size "Polynomial order must be less than the window size."
	@assert boundary_mode in (:interpolation, :nearest) "boundary_mode must be one of :interpolation, :nearest"

	# window size is 2m + 1 points
   	m = (window_size - 1) ÷ 2
	
	# build the Vandermonde design matrix A. Each row corresponds to a point in the fitting window -m:m
	# and each columns correspond to powers in the range 0:polynomial_order
   	fitting_points = -m:m
   	A = Matrix{Float64}(undef, window_size, polynomial_order + 1)
   	for i in 1:window_size, j in 1:polynomial_order + 1
        A[i,j] = fitting_points[i]^(j - 1)
    end

	if boundary_mode == :interpolation
		# for interpolation we'll want the full pseudo-inverse so we can calculate all the fit values at the edges
		# Ap = y
		C = pinv(A)

		# the filter coefficients are the rows of `C`
		filter_coeffs = C[deriv_order + 1,:] * factorial(deriv_order)

		# convolve with the filter coefficients with a couple extra steps:
		# 1. because of convolution will reverse coefficients we flip before
		# 2. c = conv(a,b) will return a vector of length(c) = length(a) + length(b) - 1 so we chop off the first and last m points
		smoothed = DSP.conv(reverse(filter_coeffs), y)[m+1:end-m]

		# for interpolation edge handling calculate the full fits
		if deriv_order == 0
			# if we are just smoothing then we can use the design and coefficient matrix as is
			AC = A*C
			smoothed[1:m] = (AC*y[1:window_size])[1:m]
			smoothed[end-m+1:end] = (AC*y[end-window_size+1:end])[end-m+1:end]
		else
			# otherwise we need to differentiate the polynomial coefficients
			# first m points
			p = C * y[1:window_size]
			for _ in 1:deriv_order
				p = [(i-1)*p[i] for i in 2:lastindex(p)]
			end
			smoothed[1:m] = A[1:m, 1:size(A,2)-deriv_order]*p
			# last m points
			p = C * y[end-window_size+1:end]
			for _ in 1:deriv_order
				p = [(i-1)*p[i] for i in 2:lastindex(p)]
			end
			smoothed[end-m+1:end] = A[m+2:end, 1:size(A,2)-deriv_order]*p

		end

		return smoothed

	elseif boundary_mode == :nearest

		# here we only need a single set of coefficients and so we can least-squares solve AᵀCᵀ = I for for a single row by picking a single column out of I
		Icol = zeros(Float64, polynomial_order+1, 1)
		Icol[deriv_order + 1] = 1.0
		filter_coeffs = transpose(A) \ Icol

		# pad the signal with the endpoints
		padded_y = [y[1] * ones(m); vec(y); y[end] * ones(m)]

		# convolve with filter
		smoothed = DSP.conv(filter_coeffs[end:-1:1], padded_y)

		# and return the valid midsection
		return smoothed[window_size:end-2*m]

	end
end


