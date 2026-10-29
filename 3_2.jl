
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
    s = s_bar * ((100.0 - A - t/365.0)/100.0)^1.8
    k_ν = k_ν_bar * (1 + t/t3)*(t+t1)/(t+t2)
    e = e_bar * (1 + t/t3)*(t+t1)/(t+t2)
    k_νm = k_νm_bar * (1 + t/t3)  # *(t+t1)/(t+t2)
    k_Ab = k_Ab_bar * (1 + t/t3)*(t+t1)/(t+t2)
    
    dT4n_spc = ν*f*(1 - e*V/(c₁+V))*s + (1 - ϕ)*ρ₄*V/(c₂+V)*T4n_spc - k_ν*V/(c₂+V)*T4n_spc - μ₄ₙ*T4n_spc
    dT4e_spc = ϕ*ρ₄*V/(c₂+V)*T4n_spc - k_ν*V/(c₂+V)*T4e_spc - μ₄ₑ*T4e_spc
    dT4n_spc_i = ν*f*e*V/(c₁+V)*s + ((1 - ϕ)*(1 - p) - p)*ρ₄*V/(c₂+V)*T4n_spc_i + k_ν*V/(c₂+V)*T4n_spc - (k₈*δ*T8e_spc + μⁱ₄ₙ)*T4n_spc_i
    dT4e_spc_i = ϕ*(1 - p)*ρ₄*V/(c₂+V)*T4n_spc_i + k_ν*V/(c₂+V)*T4e_spc - (k₈*δ*T8e_spc + μⁱ₄ₑ)*T4e_spc_i
    
    # independent from above
    dT4n_nsp = (1 - ν)*f*(1 - e*V/(c₁+V))*s + (1 - ϕ)*r₄*T4n_nsp - k_ν*V/(c₂+V)*T4n_nsp - μ₄ₙ*T4n_nsp
    dT4e_nsp = ϕ*r₄*T4n_nsp - k_ν*V/(c₂+V)*T4e_nsp - μ₄ₑ*T4e_nsp
    dT4n_nsp_i = (1 - ν)*f*e*V/(c₁+V)*s + ((1 - ϕ)*(1 - p) - p)*r₄*T4n_nsp_i + k_ν*V/(c₂+V)*T4n_nsp - (k₈*δ*T8e_spc + μⁱ₄ₙ)*T4n_nsp_i
    dT4e_nsp_i = ϕ*(1 - p)*r₄*T4n_nsp_i + k_ν*V/(c₂+V)*T4e_nsp - (k₈*δ*T8e_spc)*T4e_nsp_i

    # not independent from above, is dependent from R
    dT8n_spc = ν*(1 - f)*s + (1 - ϕ)*R*T8n_spc - μ₈ₙ*T8n_spc
    dT8e_spc = ϕ*R*T8n_spc - μ₈ₑ*T8e_spc

    # not used anywhere
    dT8n_nsp = (1 - ν)*(1 - f)*s + (1 - ϕ)*r₈*T8n*nsp - μ₈ₙ*T8n_nsp  # can be omitted
    dT8e_nsp = ϕ*r₈*T8n_nsp - μ₈ₑ*T8e_nsp  # can be omitted

    dT_ii = p*(r₄*T4n_nsp_i + ρ₄*V/(c₂+V)*T4n_spc_i) - (k₈*T8e_spc + μii)*T_ii
    dS_ab = ρAb*(T4e_spc + T4e_spc_i)*V/(c₂+V) - k_Ab*V*S_ab - μAb*S_ab
    dM_i = k_νm*V/(c₂+V) * (M - M_i) - (k₈*δ*T8e_spc + μM)*M_i
    dV = pᵢₘ*M_i/aM + N*μii*T_ii/aT - k_m * V/(c₂+V) * M/aM - e*V/(c₁+V)*s/aT - k_ν/aT*V/(c₂+V) * (T4n_spc + T4e_spc + T4n_nsp + T4e_nsp) - k_Ab*S_ab*V - μV*V 

    du .= [dT4n_spc, dT4e_spc, dT4n_spc_i, dT4e_spc_i, dT4n_nsp, dT4e_nsp, dT4n_nsp_i, dT4e_nsp_i, dT8n_spc, dT8e_spc, dT8n_nsp, dT8e_nsp, dT_ii, dS_ab, dM_i, dV]

end

#=
A, ν, f, e_bar, s_bar, ϕ, c₁, c₂, k_ν_bar = 36, 1e-3, 0.524, 0.064, 6.09, 0.64, 616.6, 1000.0, 0.065
ρ₄, ρ₈, ρAb, k₈, k_m, k_νm_bar = 1.98, 0.36, 155.0, 2.5, 60, 1.19
μ₄ₙ, μ₄ₑ, μ₈ₙ, μ₈ₑ, μii, μM, μV, μAb, k_Ab_bar, p, pᵢₘ = 0.005, 0.015, 0.006, 0.018, 0.470, 0.087, 3.0, 0.023, 0.0001, 0.03, 34.0
r₄, r₈, M, N = 0.0097, 0.0091, 360, 850
=#

c₂ = 1000.0
c₁ = 616.6
δ = 0.001
e = 0.064
f = 0.524
k₈ = 2.5
k_m = 60
k_ν_bar = 0.089
k_νm_bar = 1.19
M = 360
μAb = 0.023
μM = 0.087
μ₄ₑ = 0.015
μ₄ₙ = 0.005
μ₈ₑ = 0.018
μ₈ₙ = 0.006
μii = 0.470
μV = 3.0
N = 850
p = 0.03
A = 36
ϕ = 0.64
pᵢₘ = 34
ρAb = 155
ρ₄ = 1.98
ρ₈ = 0.36 
r₄ = 0.0097
r₈ = 0.0091
ν = 1e-3
s_bar = 6.09
τ = 7


params = A, ν, f, e_bar, s_bar, ϕ, c₁, c₂, k_ν_bar, ρ₄, ρ₈, ρAb, k₈, k_m, k_νm_bar, μ₄ₙ, μ₄ₑ, μ₈ₙ, μ₈ₑ, μii, μM, μV, μAb, k_Ab_bar, p, pᵢₘ, r₄, r₈, M, N

T4n_spc0, T4e_spc0, T4n_spc_i0, T4e_spc_i0, T4n_nsp0, T4e_nsp0, T4n_nsp_i0, T4e_nsp_i0, T_ii0, T8n_spc0, T8e_spc0, T8n_nsp0, T8e_nsp0, S_ab0, M_i0, V0 = 0.92, 0.0, 0.0, 0.0, 920.0, 379.0, 0.0, 0.0, 0.0, 0.47, 0.0, 467.0, 150.0, 0.0, 0.0, 3.0

