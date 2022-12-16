using DifferentialEquations, CairoMakie, UnPack, DataDrivenDiffEq, StaticArrays

struct parameter_stauffer
    Ω :: Float64 
    φ :: Float64
    A :: Float64
    k :: Float64
    Q :: Float64
    d :: Float64 
    C₂ :: Float64
    C₈ :: Float64 
    function parameter_stauffer(amplitude, quality, distance)
        new(1., 0., amplitude, 1., quality, distance, 0.3, 1.)
    end
end

struct parameter_lee
    A₁ :: Float64
    A₂ :: Float64
    R :: Float64
    Q :: Float64
    E :: Float64
    k :: Float64
    A :: Float64
    f₀ :: Float64
    φ :: Float64
    d :: Float64
    Ω :: Float64
    function parameter_lee(amplitude, f₀, φ, distance) 
        new(1.35961e-70, 1.8651e-19, 150e-9, 100., 176e9, 40., amplitude, f₀, φ, distance, 1.)
    end
end

function lennard_jones_lee(u, p, t)
    @unpack A₁, A₂, R, Q, E, k, A, f₀, φ, d, Ω = p
    z = (u[1] - d)
    dx = u[2]
    dy = -1/Q*u[2] - k*(u[1] - d) - A₁*R/(180. * z^8) + A₂*R/(6. * z^2) + A*sin(Ω*t + φ)
    SA[dx, dy] 
end

function lennard_jones_stauffer(u, p, t)
    @unpack Ω, φ, A, k, Q, d, C₂, C₈ = p
    dx = u[2] 
    dy = -1/Q * u[2] - (u[1] - d) + A*sin(Ω*t + φ) - C₈ * u[1]^(-8)+ C₂ * u[1]^(-2)
    SA[dx, dy]
end


p = parameter_stauffer(80e-9, 500., 100e-9)
p_lee = parameter_lee(1e-9, 3., 0., 100e-9)
Δt = 0.005
tspan = (0., 4000.)
u0 = SA[0., 0.]
prob = ODEProblem(lennard_jones_stauffer, u0, tspan, p)
prob_lee = ODEProblem(lennard_jones_lee, u0, tspan, p_lee)
sol_lee = solve(prob_lee, AutoTsit5(Rosenbrock23()), dt=Δt, adaptive=false)
