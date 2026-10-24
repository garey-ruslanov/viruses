
using DifferentialEquations, Plots

function rightside!(du, u, P, t)
    du .= P.A * u + P.b
end

u0 = [10.0, 5.0, 0.0]
P = (A = [-2.0 1.0 0.0; 1.0 -2.0 1.0; 0.0 1.0 -2.0], b = [-0.5, -0.5, -0.5])
timespan = (0.0, 10.0)
problem = ODEProblem(rightside!, u0, timespan, P)
