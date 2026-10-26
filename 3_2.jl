
using DifferentialEquations, Plots

function rightside!(du, u, params, t)
    A, ν, f, e_bar, s_bar, ϕ, c₁, c₂, k_ν_bar, ρ₄, ρ₈, ρAb, k₈, k_m, k_νm_bar, μ₄ₙ, μ₄ₑ, μ₈ₙ, μ₈ₑ, μii, μM, μV, μAb, k_Ab_bar, p, pᵢₘ, r₄, r₈, M, N = params

    @views T4n_spc, T4e_spc, T4n_spc_i, T4e_spc_i, T4n_nsp, T4e_nsp, T4n_nsp_i, T4e_nsp_i, T8n_spc, T8e_spc, T8n_nsp, T8e_nsp, T_ii, S_ab, M_i, V = u

    t1, t2, t3 = 7, 42, 40000
    aT = 0.06
    aM = 0.03
    μⁱ₄ₑ = μ₄ₑ
    μⁱ₄ₙ = μ₄ₙ
    δ = 0.01

    R = ρ₈*(T_ii + δ*(T4n_nsp_i + T4n_spc_i + T4e_nsp_i + T4e_spc_i + aT/aM * M_i))

    

end

#=
A, ν, f, e_bar, s_bar, ϕ, c₁, c₂, k_ν_bar = 36, 1e-3, 0.524, 0.064, 6.09, 0.64, 616.6, 1000.0, 0.065
ρ₄, ρ₈, ρAb, k₈, k_m, k_νm_bar = 1.98, 0.36, 155.0, 2.5, 60, 1.19
μ₄ₙ, μ₄ₑ, μ₈ₙ, μ₈ₑ, μii, μM, μV, μAb, k_Ab_bar, p, pᵢₘ = 0.005, 0.015, 0.006, 0.018, 0.470, 0.087, 3.0, 0.023, 0.0001, 0.03, 34.0
r₄, r₈, M, N = 0.0097, 0.0091, 360, 850
=#

params = A, ν, f, e_bar, s_bar, ϕ, c₁, c₂, k_ν_bar, 
          ρ₄, ρ₈, ρAb, k₈, k_m, k_νm_bar, 
          μ₄ₙ, μ₄ₑ, μ₈ₙ, μ₈ₑ, μii, μM, μV, μAb, k_Ab_bar, 
          p, pᵢₘ, r₄, r₈, M, N
