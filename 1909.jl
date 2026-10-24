
using DifferentialEquations, Plots

function h_constant(p, t; idxs=1)
    t ≤ 0 || error("fuck off")
    return 0.0
end

function dde_constant_delay_trivial!(du, u, h, p, t)
    τ = p
    du[1] = -1*h(p, t - τ; idxs=1)
    return nothing
end

function dde_logistical!(du, u, h, p, t)
    τ, b, C = p
    du[1] = b*u[1]*(1 - h(p, t - τ, idxs=1)/C)

    return nothing
end

τ_array = [0.2, 1.4, 1.8, 3.6]
τ = 0.2
b = 1.0
C = 1.0
u0 = [0.2]

let plt = plot()
for τ in τ_array
    p = τ, b, C
    timespan = (0.0, 60.0)
    problem = DDEProblem(dde_logistical!, u0, h_constant, timespan, p; constant_lags=[τ])
    solution = solve(problem, MethodOfSteps(Tsit5()); abstol=1e-10, reltol=1e-8)
    plot!(plt, solution; label="τ=$τ")
end
display(plt)
savefig(plt, "12314560.png")
end