using CairoMakie, LinearAlgebra

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

Ω = 1.5
x = collect(0.:0.005:2π)
y = sin.(Ω*x .+ π/5) .+ cos.(2*Ω*x)

base = sinusoidal_bases(3, Ω, 0., 2π)
θ = lsq_regression(x, y, base)
ϕ_s = atan(θ[2]/θ[3])
ϕ_c = atan(θ[4]/θ[5])

fig = Figure(resolution=(800, 800))
ax = Axis(fig[1, 1])
lines!(ax, x, y)
lines!(ax, x, cos.(Ω*x .- ϕ_s) .+ cos.(2*Ω*x .- ϕ_c))
fig

