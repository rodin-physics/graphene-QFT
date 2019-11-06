using QuadGK
using Plots
using LinearAlgebra
using LaTeXStrings
using DelimitedFiles

## Parameters
const ν = 1e-6;       # Small number for relative tolerance
const η = 1e-8;       # Small number for absolute tolerance
const NumEvals = 1e6; # Maximum number of evaluations in integrals

# Graphene hopping integral and lattice vectors in Angstroms
const t = 2.8;
const d1 = 2.46 .* [1 √(3)] ./ 2;
const d2 = 2.46 .* [-1 √(3)] ./ 2;

# Sublattice spinors
const A = [1, 0]
const B = [0, 1]

# Graphene Coordinate type: R = u * d1 + v * d2
struct GrapheneCoord
    u :: Int
    v :: Int
    sublattice :: Vector{Int64}
end

## Helper functions
# Auxiliary W function
function W(z :: ComplexF64, x :: Float64)
    return (((z / t)^2 - 1.0) / (4.0 * cos(x)) - cos(x))
end

# Auxiliary Y function used in the calculation of Ω
function Y(n :: Int, z, x)
    W_ = W(z, x)
    return (W_ - √(W_ - 1) * √(W_ + 1))^n / (√(W_ - 1) * √(W_ + 1))
end

# When computing Ω, occasionally the integrand becomes small enough to give NaN
# This helper functions is used to catch these instances
function Ω_Integrand(z, u, v, x :: Float64)
    res = exp(1.0im * (u - v) * x) / cos(x) * Y(abs.(u + v), z, x)
    if isnan(res) == true
        return 0.0 + 0.0im
    else
        return res
    end
end

## Main functions

function Ω(z, u, v)
    f_int(x) = Ω_Integrand(z, u, v, x)
    res = quadgk(f_int, 0, 2 * π, maxevals = NumEvals, atol = η)
    return (res[1] / (8.0 * π * t^2) )
end

# Atom self-energy function
function Σ(z)
    return (z * Ω(z, 0, 0))
end

# The propagator function picks out the correct element of the Ξ matrix based
# on the sublattices of the graphene coordinates
function Propagator(Imp_l :: GrapheneCoord, Imp_m :: GrapheneCoord, z)
    u = Imp_l.u - Imp_m.u
    v = Imp_l.v - Imp_m.v
    if Imp_l.sublattice == Imp_m.sublattice
        return (z * Ω(z, u, v))
    elseif ([Imp_l.sublattice, Imp_m.sublattice] == [A, B])
        return (- t * (Ω(z, u, v) + Ω(z, u + 1, v) + Ω(z, u, v + 1)))
    elseif ([Imp_l.sublattice, Imp_m.sublattice] == [B, A])
        return (- t * (Ω(z, u, v) + Ω(z, u - 1, v) + Ω(z, u, v - 1)))
    else
        error("Illegal sublattice parameter")
    end
end

# The (I^T Ξ I) Matrix
function Propagator_Matrix(z, Coords :: Vector{GrapheneCoord})
    nCoords = length(Coords)

    CoordsT_Mat = repeat(Coords, 1, nCoords)  # Position matrix
    Coords_Mat = permutedims(CoordsT_Mat)     # Transpose of position matrix

    Prop = (map((x, y) -> Propagator(x, y, z), CoordsT_Mat, Coords_Mat))
    return Prop
end


## Data Processing Functions

# In order to plot the calculated results as a lattice, one needs to create
# the correct coordinate arrays

function Data_Process(A_Lattice, B_Lattice)
    sz = size(A_Lattice)
    nPts = floor(Int,(sz[1] - 1)/2);
    # Lattice shift between the two sublattices
    latticeShift = -1/√(3) * 2.46;
    # Arrays
    d1s = -nPts : 1 : nPts; # d1 vectors
    d2s = -nPts : 1 : nPts; # d2 vectors

    D1S = repeat(d1s, 1, 2 * nPts + 1);
    D2S = repeat(d2s', 2 * nPts + 1, 1);

    # Coordinates of the carbon atoms
    XS = d1[1] .* D1S .+ d2[1] .* D2S;
    YS = d1[2] .* D1S .+ d2[2] .* D2S;

    # Flatten the coordinates and the data
    XS_A = reshape(XS, 1, sz[1]^2);
    YS_A = reshape(YS, 1, sz[1]^2);
    A_Lattice = reshape(A_Lattice, 1, sz[1]^2);

    XS_B = reshape(XS, 1, sz[1]^2);
    YS_B = reshape(YS, 1, sz[1]^2) .+ latticeShift;
    B_Lattice = reshape(B_Lattice, 1, sz[1]^2);

    # Combine all the coordinates and data for both sublattices
    XS = vcat(XS_A', XS_B');
    YS = vcat(YS_A', YS_B');
    dta = vcat(A_Lattice', B_Lattice');

    return([XS YS dta])
end
