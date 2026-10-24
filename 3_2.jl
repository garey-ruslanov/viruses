
using DifferentialEquations, Plots

function rightside!(du, u, params, t)
    A, ν, f, e_bar, s_bar, ϕ, c₁, c₂, k_ν_bar, ρ₄, ρ₈, ρAb, k₈, k_m, k_νm_bar, μ₄ₙ, μ₄ₑ, μ₈ₙ, μ₈ₑ, μii, μM, μV, μAb, k_Ab_bar, p, pᵢₘ, r₄, r₈, M, N = params

    

    t1, t2, t3 = 7, 42, 40000
    aT = 0.06  # 0.06
    aM = 0.03  # 0.03
    μⁱ₄ₑ = μ₄ₑ
    μⁱ₄ₙ = μ₄ₙ
    δ = 0.01

    

end

A, ν, f, e_bar, s_bar, ϕ, c₁, c₂, k_ν_bar = 36, 1e-3, 0.524, 0.064, 6.09, 0.64, 616.6, 1000.0, 0.065
ρ₄, ρ₈, ρAb, k₈, k_m, k_νm_bar = 1.98, 0.36, 155.0, 2.5, 60, 1.19
μ₄ₙ, μ₄ₑ, μ₈ₙ, μ₈ₑ, μii, μM, μV, μAb, k_Ab_bar, p, pᵢₘ = 0.005, 0.015, 0.006, 0.018, 0.470, 0.087, 3.0, 0.023, 0.0001, 0.03, 34.0
r₄, r₈, M, N = 0.0097, 0.0091, 360, 850

params = A, ν, f, e_bar, s_bar, ϕ, c₁, c₂, k_ν_bar, 
          ρ₄, ρ₈, ρAb, k₈, k_m, k_νm_bar, 
          μ₄ₙ, μ₄ₑ, μ₈ₙ, μ₈ₑ, μii, μM, μV, μAb, k_Ab_bar, 
          p, pᵢₘ, r₄, r₈, M, N
