using LinearAlgebra
"""
This algorithm is used to get information of the phase and amplitude at each
timestep. --> extremely fast simulations
"""

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

function update_basis!(b::AbstractArray{Float64}, bases::AbstractArray{Function}, t::Float64)
    b .= [f(t) for f in bases]
end


function update_weights!(w::AbstractArray{Float64}, b::AbstractArray{Float64}, data::Float64, μ::Float64)
    w .= w .+ μ .* inv(dot(b, b)) .* b .* ϵ_lms(data, w, b)
end

b = zeros(7)
w = zeros(7)
base = create_bases(3, 1)
Δt = 0.005
ts =collect(0.:Δt:30π)
y = 2*sin.(ts) + 3*cos.(2*ts)

for idx in eachindex(y)
    update_basis!(b, base, ts[idx])
    update_weights!(w, b, y[idx], Δt)
end

# works wonderful, needs some periods to converge