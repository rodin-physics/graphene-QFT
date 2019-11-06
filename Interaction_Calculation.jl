using Distributed
using DelimitedFiles

@everywhere include("Computation_Settings.jl")

# Calculation
@everywhere function f_A_Sublattice(x1, x2)
    println("x1 = $(x1), x2 = $(x2)")
    Loc = GrapheneCoord(x1, x2, A)
    res = F_Dopants(μ, [Loc_Origin, Loc], τ, ϵ)
    println(res)
    return res[1]
end

@everywhere function f_B_Sublattice(x1, x2)
    println("x1 = $(x1), x2 = $(x2)")
    Loc = GrapheneCoord(x1, x2, B)
    res = F_Dopants(μ, [Loc_Origin, Loc], τ, ϵ)
    println(res)
    return res[1]
end

resA = pmap(f_A_Sublattice, D1S, D2S)
resB = pmap(f_B_Sublattice, D1S, D2S)

singleAdsorbate = F_Dopants(μ, [Loc_Origin], τ, ϵ)
resA = resA .- 2 * singleAdsorbate
resB = resB .- 2 * singleAdsorbate

writedlm("Data/Interaction/F_A_mu_02_A.txt", resA)
writedlm("Data/Interaction/F_A_mu_02_B.txt", resB)
