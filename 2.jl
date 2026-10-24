
# 2 method

function rightside!(du, u, p, t)
    τPimm, τPmat, τRimm, τRmat, τC, IP, IC, Pimm0, Pmat0, Rimm0, Rmat0, a0, f, α, ε, γ, θ, ν, h, L, λ, cPmat, ρR, ρC, ζ = p
    #λ = 1.0
    #cPmat = 20.0
    #ρR = 0.0
    #ρC = 0.0
    #ζ = 0.0
    @views Pimm, Pmat, Rimm, Rmat, a, C = u
    function constrain(at, L, a0, h)
        if log(at / a0) < L
            return 1.0
        else
            return h*log(at / a0)
        end
    end
    dPimm = (IP + f*(Pmat0 - Pmat + Rmat0 - Rmat)) / constrain(a, L, a0, h) - τPimm*Pimm - 0.0
    dPmat = τPimm*Pimm - τPmat*Pmat - cPmat*a*C*Pmat
    dRimm = 2/3* (IP + f*(Pmat0 - Pmat + Rmat0 - Rmat)) / constrain(a, L, a0, h) - τRimm*Rimm
    dRmat = τRimm*Rimm - (τRmat - ρR)*Rmat
    da = a*(θ - ζ - γ*C)
    dC = λ*a*(ε*IC + α*C)*(Pmat/Pmat0)^ν - (τC - ρC)*C
    du .= [dPimm, dPmat, dRimm, dRmat, da, dC]
end

Pimm0, Pmat0, Rimm0, Rmat0, a0, C0 = 5.0, 100.0, 3.33, 66.7, 5e-4, 0.0
τPimm, τPmat, τRimm, τRmat, τC, IP, IC = 0.2, 0.01, 0.2, 0.01, 0.01, 1.0, 0.2
f, α, ε, γ, θ, ν, h, L = 0.01, 0.7, 0.512, 0.3, 0.02, 1.6, 0.3, 3.0
λ, cPmat, ρR, ρC, ζ = 1.0, 20.0, 0.0, 0.0, 0.005
params = τPimm, τPmat, τRimm, τRmat, τC, IP, IC, Pimm0, Pmat0, Rimm0, Rmat0, a0, f, α, ε, γ, θ, ν, h, L, λ, cPmat, ρR, ρC, ζ

u0 = [Pimm0, Pmat0, Rimm0, Rmat0, a0, C0]
timespan = (0, 1500)
problem = ODEProblem(rightside!, u0, timespan, params)
solution = solve(problem, Rosenbrock23())

plo = plot(solution, idxs=(0, 2), label="CD4 mature")
plot!(plo, solution, idxs=(0, 4), label="CD8 mature")
savefig(plo, "cd4cd8mature_zeta_5e_3.png")
plo = plot(solution, idxs=(0, 1), label="CD4 immature")
plot!(plo, solution, idxs=(0, 3), label="CD8 immature")
savefig(plo, "cd4cd8immature_zeta_5e_3.png")

plo = plot(solution, idxs=(0, 6), label="Cytotoxic CD8")
savefig(plo, "ctoxcd8_zeta_5e_3.png")

plo = plot(solution, idxs=(0, 5), label="Virus RNA", yscale=:log10)
savefig(plo, "virus_zeta_5e_3.png")