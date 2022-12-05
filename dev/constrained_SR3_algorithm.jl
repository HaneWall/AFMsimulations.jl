using LinearAlgebra 

"""
Problem we would like to optimize: 

min_{Ξ, W} 1/2 ||Θ(X)Ξ - d/dtX||^2 + R(W, λ, ν) + 1/(2ν) ||Ξ - W||^2 . 

sometimes even more complicated, we would like to solve this problem with
constraints. We write the constraints in the matrix form C Ξ = d. 

X .. Data Points of the position at different timeshots 
Θ .. library/candidate fucntions
Ξ .. full coefficients 
W .. sparse coefficients 
R .. thresholding operator

Let us first seek the solution for no constraints. Here we can use the normal
equation form of the problem and define variables to seek a least-square
solution. Remember: min_{x} ||Ax - b||^2 iff A'Ax = A'b. A'A can then be
factored by a cholesky decomposition. 

---> 1. Least squares

H = (Θ'Θ + 1/ν I)^{-1}
H = ΘΘ' + νI [ == (A'A)^{-1}
X = cholesky(H) --> upper triangular form of H 
Y = dX/dt Θ' [ == bA']
Ξ = Y / X   
=> Ξ * (Θ + 1/ν I) ≈ dX/dt  


---> 2. Apply Regularization / Sparsity Promotion  
Now that we got the solution for the least squares problem, we need to
update W and R(W). 

W = prox(Ξ) --> sparsify
active_set! --> below threshhold candidate is no longer in computational domain. 
active_set! applies a threshhold via a chosen proximity operator (hard threshholding in our case)

---> 3. define proximity operator
There exist many different proximity operator, Stauffer et al use the l_0 norm
which corresponds to hard thresholding. λ = η^2/(2ν), where η is the desired
hard threshhold. As mentioned in Stauffer we can use η=1/Q as a good starting point.

"""

"""
returns still active candidate functions (to be exact their boolean index) after each regression
"""
function active_set!(idx::BitArray, x::AbstractArray{T}, λ::T) where {T}
@assert size(idx) == size(x)
@inbounds foreach(eachindex(x)) do i
    idx[i] = abs(x[i]) > sqrt(2 * λ)
    end
    return
end