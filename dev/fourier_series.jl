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
y = 1/sqrt(2) .* cos.(Ω*x .+ π/6) .+  1.5 .* sin.(2*Ω * x)

base = sinusoidal_bases(3, Ω, 0., 2π)
θ = lsq_regression(x, y, base)
A_s = sqrt(θ[2]^2 + θ[3]^2) 
ϕ_s = atan(θ[3], θ[2])

A_c = sqrt(θ[4]^2 + θ[5]^2) 
ϕ_c = atan(θ[5], θ[4])

fig = Figure(resolution=(800, 800))
ax = Axis(fig[1, 1])
lines!(ax, x, y, linewidth=3.5)
lines!(ax, x, A_s .* sin.(Ω*x .+ ϕ_s) .+ A_c .*sin.(2*Ω*x .+ ϕ_c) , linewidth=3.5, linestyle=:dash)
fig

