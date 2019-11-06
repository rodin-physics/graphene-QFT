using Profile
include("Lib/general.jl")
include("Lib/graphene_dopants.jl")

# Parameters

ϵ = 0.35;       # Energy of the host atom
τ = t;          # Change in coupling between the host atoms and their neighbors
                # NOTE: Pristine graphene has hopping equal to -t
μ = -0.2;       # Chemical potential

nPts = 25;     # Number of grid points away from zero

# Arrays
d1s = -nPts : 1 : nPts; # d1 vectors
d2s = -nPts : 1 : nPts; # d2 vectors

# Lattice vector matrices
D1S = repeat(d1s, 1, 2 * nPts + 1);
D2S = repeat(d2s', 2 * nPts + 1, 1);

# Impurity Arrangement
Loc_Origin = GrapheneCoord(0, 0, A)
Loc_1 = GrapheneCoord(10, -10, A)
Loc_2 = GrapheneCoord(-10, 10, A)

Dopants = [Loc_1, Loc_2]
