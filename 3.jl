
# 3 method

function rightside!(du, u, params, t)
    A, ν, f, e_bar, s_bar, ϕ, c₁, c₂, k_ν_bar, ρ₄, ρ₈, ρAb, k₈, k_m, k_νm_bar, μ₄ₙ, μ₄ₑ, μ₈ₙ, μ₈ₑ, μii, μM, μV, μAb, k_Ab_bar, p, pᵢₘ, r₄, r₈, M, N = params

    T4spec_n, T4spec_e, T4spec_ni, T4spec_ei, T4nspc_n, T4nspc_e, T4nspc_ni, T4nspc_ei, T8spec_n, T8spec_e, T8nspc_n, T8nspc_e, Tii, SAb, Mi, V = u
    
    t1, t2, t3 = 7, 42, 40000
    aT = 0.06  # 0.06
    aM = 0.03  # 0.03
    μⁱ₄ₑ = μ₄ₑ
    μⁱ₄ₙ = μ₄ₙ
    δ = 0.01
    
    ut = (t + t1)/(t + t2)
    R = ρ₈*(Tii + δ*(T4spec_ni + T4nspc_ni + T4spec_ei + T4nspc_ei + aT/aM*Mi))
    s = s_bar*((100 - A - t/365) / 100)^1.8
    k_ν = k_ν_bar*(1 + t/t3)*ut
    e = min(1.0, e_bar*(1 + t/t3)*ut)
    k_νm = k_νm_bar*(1 + t/t3)*ut
    k_Ab = k_Ab_bar*(1 + t/t3)*ut

    dT4spec_n = ν*f*(1 - e*V/(c₁+V))*s + (1 - ϕ)*(ρ₄*V)/(c₂+V)*T4spec_n - k_ν*V/(c₂+V)*T4spec_n - μ₄ₙ*T4spec_n
    dT4spec_e = ϕ*ρ₄*V/(c₂+V)*T4spec_n - k_ν*V/(c₂+V)*T4spec_n - μ₄ₑ*T4spec_e
    dT4spec_ni = ν*f*(e*V)/(c₁+V)*s + ((1-ϕ)*(1-p) - p)*ρ₄*V/(c₂+V)*T4spec_ni + k_ν*V/(c₂+V)*T4spec_n - (k₈*δ*T8spec_e + μⁱ₄ₙ)*T4spec_ni
    dT4spec_ei = ϕ*(1-p)*ρ₄*V/(c₂+V)*T4spec_ni + k_ν*V/(c₂+V)*T4spec_e - (k₈*δ*T8spec_e + μⁱ₄ₑ)*T4spec_ei

    dT4nspc_n = (1 - ν)*f*(1 - e*V/(c₁+V))*s + (1 - ϕ)*r₄*T4nspc_n - k_ν*V/(c₂+V)*T4nspc_n - μ₄ₙ*T4nspc_n
    dT4nspc_e = ϕ*r₄*T4nspc_n - k_ν*V/(c₂+V)*T4nspc_e - μ₄ₑ*T4nspc_e
    dT4nspc_ni = (1 - ν)*f*e*V/(c₁+V)*s + ((1-ϕ)*(1-p) - p)*r₄*T4nspc_ni + k_ν*V/(c₂+V)*T4nspc_n - (k₈*δ*T8spec_e + μⁱ₄ₙ)*T4nspc_ni
    dT4nspc_ei = ϕ*(1-p)*r₄*T4nspc_ni + k_ν*V/(c₂+V)*T4nspc_e - (k₈*δ*T8spec_e + μⁱ₄ₑ)*T4nspc_ei

    dT8spec_n = ν*(1 - f)*s + (1 - ϕ)*R*T8spec_n - μ₈ₙ*T8spec_n
    dT8spec_e = ϕ*R*T8spec_n - μ₈ₑ*T8spec_e

    dT8nspc_n = (1 - ν)*(1 - f)*s + (1 - ϕ)*r₈*T8nspc_n - μ₈ₙ*T8nspc_n
    dT8nspc_e = ϕ*r₈*T8nspc_n - μ₈ₑ*T8nspc_e

    dTii = p*(r₄*T4nspc_ni + ρ₄*V/(c₂+V)*T4spec_ni) - (k₈*T8spec_e + μii)*Tii
    dSAb = ρAb*(T4spec_e + T4spec_ei)*V/(c₂+V) - k_Ab*V*SAb - μAb*SAb
    dMi = k_νm*V/(c₂+V)*(M - Mi) - (k₈*δ*T8spec_e + μM)*Mi
    dV = pᵢₘ*Mi/aM + N*μii*Tii/aT - k_m*V/(c₂+V)*M/aM - e*V/(c₁+V)*s/aT - k_ν/aT*V/(c₂+V)*(T4spec_n + T4spec_e + T4nspc_n + T4nspc_e) - k_Ab*SAb*V - μV*V

    du .= [dT4spec_n, dT4spec_e, dT4spec_ni, dT4spec_ei, dT4nspc_n, dT4nspc_e, dT4nspc_ni, dT4nspc_ei, dT8spec_n, dT8spec_e, dT8nspc_n, dT8nspc_e, dTii, dSAb, dMi, dV]
end

A, ν, f, e_bar, s_bar, ϕ, c₁, c₂, k_ν_bar = 36, 1e-3, 0.524, 0.064, 6.09, 0.64, 616.6, 1000.0, 0.065
ρ₄, ρ₈, ρAb, k₈, k_m, k_νm_bar = 1.98, 0.36, 155.0, 2.5, 60, 1.19
μ₄ₙ, μ₄ₑ, μ₈ₙ, μ₈ₑ, μii, μM, μV, μAb, k_Ab_bar, p, pᵢₘ = 0.005, 0.015, 0.006, 0.018, 0.470, 0.087, 3.0, 0.023, 0.0001, 0.03, 34.0
r₄, r₈, M, N = 0.0097, 0.0091, 360, 850

params = A, ν, f, e_bar, s_bar, ϕ, c₁, c₂, k_ν_bar, 
          ρ₄, ρ₈, ρAb, k₈, k_m, k_νm_bar, 
          μ₄ₙ, μ₄ₑ, μ₈ₙ, μ₈ₑ, μii, μM, μV, μAb, k_Ab_bar, 
          p, pᵢₘ, r₄, r₈, M, N

T4spec_n0, T4spec_e0, T4spec_ni0, T4spec_ei0, T4nspc_n0, T4nspc_e0, T4nspc_ni0, T4nspc_ei0, T8spec_n0, T8spec_e0, T8nspc_n0, T8nspc_e0, Tii0, SAb0, Mi0, V0 = 0.92, 0.0, 0.0, 0.0, 920.0, 379.0, 0.0, 0.0, 0.47, 0.0, 467.0, 150.0, 0.0, 0.0, 0.0000001, 3.0

u0 = [T4spec_n0, T4spec_e0, T4spec_ni0, T4spec_ei0, T4nspc_n0, T4nspc_e0, T4nspc_ni0, T4nspc_ei0, T8spec_n0, T8spec_e0, T8nspc_n0, T8nspc_e0, Tii0, SAb0, Mi0, V0]
timespan = (0, 1000)
problem = ODEProblem(rightside!, u0, timespan, params)
solution = solve(problem, Rodas5(), reltol=1e-8, abstol=1e-10)

#=

plo = plot(solution, idxs=(0, 9), label="T8 spec., not eff.")
savefig(plo, "T8spec_noteff_.png")
plo = plot(solution, idxs=(0, 10), label="T8 spec., eff.")
savefig(plo, "T8spec_eff_.png")
plo = plot(solution, idxs=(0, 11), label="T8 not spec., not eff.")
savefig(plo, "T8nspec_noteff_.png")
plo = plot(solution, idxs=(0, 12), label="T8 not spec., eff.")
savefig(plo, "T8nspec_eff_.png")

mask = solution[13,:] .>= 1e-8
uvals = solution[13,:][mask]
tvals = solution.t[mask]
plo = plot(tvals, uvals, yscale=:log10, label="T active infected")
savefig(plo, "T_activeinfected_log.png")

mask = solution[15,:] .>= 1e-1
uvals = solution[15,:][mask]
tvals = solution.t[mask]
plo = plot(tvals, uvals, yscale=:log10, label="Macrophage infected")
savefig(plo, "Macrophage_infected_log.png")

mask = solution[14,:] .>= 1e-1
uvals = solution[14,:][mask]
tvals = solution.t[mask]
plo = plot(tvals, uvals, yscale=:log10, label="Antibody")
savefig(plo, "Antibody_log.png")

mask = solution[16,:] .>= 1e-1
uvals = solution[16,:][mask]
tvals = solution.t[mask]
plo = plot(tvals, uvals, yscale=:log10, label="Virus dynamics")
savefig(plo, "Virus_log.png")

=#

plo = plot(solution, idxs=(0, ))