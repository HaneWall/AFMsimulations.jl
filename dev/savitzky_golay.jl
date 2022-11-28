using DSP: conv 
using LinearAlgebra: pinv
using CairoMakie
CairoMakie.activate!(type="svg")


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
		smoothed = conv(reverse(filter_coeffs), y)[m+1:end-m]

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
				p = [(i-1)*p[i] for i in 2:length(p)]
			end
			smoothed[1:m] = A[1:m, 1:size(A,2)-deriv_order]*p
			# last m points
			p = C * y[end-window_size+1:end]
			for _ in 1:deriv_order
				p = [(i-1)*p[i] for i in 2:length(p)]
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
		smoothed = conv(filter_coeffs[end:-1:1], padded_y)

		# and return the valid midsection
		return smoothed[window_size:end-2*m]

	end
end


"""
Noisy sine data for test
"""

δx = 0.1
x_set = collect(0.:δx:4π)
y_set = sin.(x_set)
noise = 0.1*randn(length(x_set))
window_size = 31
y_noised = y_set .+ noise
y_filtered = savitzky_golay_filter(y_noised, window_size, 3; deriv_order=0, boundary_mode = :interpolation)
fig = Figure(resolution=(800, 800))
ax = Axis(fig[1, 1])
lines!(ax, x_set, y_noised, color=:grey)
lines!(ax, x_set, y_filtered, linewidth=2)
lines!(ax, x_set, y_set, linewidth=2)
lines!(ax, x_set, 1/δx * savitzky_golay_filter(y_noised, window_size, 3; deriv_order=1, boundary_mode = :interpolation))
fig
