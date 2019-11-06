using Distributed

@everywhere include("Lib/graphene_dopants.jl")
@everywhere include("Computation_Settings.jl")

# Calculation
@everywhere function f_A_Sublattice(x1, x2)
    println("x1 = $(x1), x2 = $(x2)")
    Loc = GrapheneCoord(x1, x2, A)
    res = Δρ_Dopants(Loc, μ, Dopants, τ, ϵ)
    println(res)
    return res[1]
end

@everywhere function f_B_Sublattice(x1, x2)
    println("x1 = $(x1), x2 = $(x2)")
    Loc = GrapheneCoord(x1, x2, B)
    res = Δρ_Dopants(Loc, μ, Dopants, τ, ϵ)
    println(res)
    return res[1]
end

resA = pmap(f_A_Sublattice, D1S, D2S)
resB = pmap(f_B_Sublattice, D1S, D2S)

writedlm("Data/Density/rho_mu_04_Dop2_A.txt", resA)
writedlm("Data/Density/rho_mu_04_Dop2_B.txt", resB)
