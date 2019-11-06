include("general.jl")

# Function for determining the neighboring coordinates of a graphene site
function neighbors(loc :: GrapheneCoord)
    u = loc.u
    v = loc.v
    sublattice = loc.sublattice
    if sublattice == A
        return ([GrapheneCoord(u, v, B);
                 GrapheneCoord(u + 1, v, B);
                 GrapheneCoord(u, v + 1, B)])
    elseif sublattice == B
        return ([GrapheneCoord(u, v, A);
                 GrapheneCoord(u - 1, v, A);
                 GrapheneCoord(u, v - 1, A)])
    else
        error("Illegal sublattice parameter")
    end
end

# Function to get the coordinates of the dopants and their unique neighbors
function scattering_atoms(dopants :: Vector{GrapheneCoord})
    # Obtain the coordinates of all the dopants' neighbors
    all_neighbors = collect(Iterators.flatten(map(neighbors, dopants)))
    # Combine all the atoms (including the dopants) into a single array
    all_atoms = unique(vcat(dopants, all_neighbors))
    return all_atoms
end

# Function to check whether there is a modified bond between two coordinates
function bondcheck(c1 :: GrapheneCoord,
                   c2 :: GrapheneCoord,
                   dopants :: Vector{GrapheneCoord})
    return (c1 in neighbors(c2) && (c1 in dopants || c2 in dopants))
end

# Integrand used to calculate the local density
function Δρ_Dopants_Integrand(loc :: GrapheneCoord, z, dopants, τ, ϵ)
    atoms = scattering_atoms(dopants)
    n_atoms = length(atoms)
    # Arrange the atoms into a matrix
    atoms_M = repeat(atoms, 1, n_atoms)
    atoms_M_T = permutedims(atoms_M)
    # Check if there are modified bonds between pairs of atoms
    Δ_ = map((x, y) -> τ * bondcheck(x, y, dopants), atoms_M, atoms_M_T)
    # Add the on-site potential for the dopants
    Δ_ = Δ_ .+ ϵ .* Diagonal(map(x -> x in dopants, atoms))
    # Get the propagator matrix I^T Ξ I and the scattering matrix
    # (1 - I^T Ξ I Δ)^(-1)
    prop_mat = Propagator_Matrix(z, atoms)
    scattering_matrix = inv(Matrix{Int}(I, n_atoms, n_atoms) .- prop_mat * Δ_)

    PropVectorL = reshape(map(x -> Propagator(loc, x , z), atoms), (1, n_atoms))
    PropVectorR = reshape(map(x -> Propagator(x, loc , z), atoms), (n_atoms, 1))

    return real((PropVectorL * Δ_ * scattering_matrix * PropVectorR)[1])
end

# Local density function
function Δρ_Dopants(loc, μ, dopants, τ, ϵ)
    f_int(x) = Δρ_Dopants_Integrand(loc, μ + 1im * x, dopants, τ, ϵ)
    res = quadgk(f_int, 0, Inf, rtol = 1e-2)
    return (2 * res[1] / (2 * π))
end

# Integrand used to compute the interaction energy. We catch NaN and return 0
function F_Dopants_Integrand(z,
                             dopants :: Vector{GrapheneCoord},
                             τ :: Float64,
                             ϵ :: Float64)
    atoms = scattering_atoms(dopants)
    n_atoms = length(atoms)
    # Arrange the atoms into a matrix
    atoms_M = repeat(atoms, 1, n_atoms)
    atoms_M_T = permutedims(atoms_M)
    # Check if there are modified bonds between pairs of atoms
    Δ_ = map((x, y) -> τ * bondcheck(x, y, dopants), atoms_M, atoms_M_T)
    # Add the on-site potential for the dopants
    Δ_ = Δ_ .+ ϵ .* Diagonal(map(x -> x in dopants, atoms))
    # Get the propagator matrix I^T Ξ I and the scattering matrix
    # (1 - I^T Ξ I Δ)^(-1)
    prop_mat =  Propagator_Matrix(z, atoms)
    log_kernel = det(Matrix{Int}(I, n_atoms, n_atoms) .- prop_mat * Δ_)

    res = -log(log_kernel) / (2 * π)

    if isnan(res) == true
        return 0.0 + 0.0im
    else
        return res
    end
end

# Impurity interaction energy
function F_Dopants(μ, dopants, τ, ϵ)
    f_int(x) = real(F_Dopants_Integrand(μ + 1im * x, dopants, τ, ϵ))
    res = quadgk(f_int, 0, Inf, rtol = ν)
    return (2 * res[1])
end

# Local spectral function
function spectralFn(loc :: GrapheneCoord,
                    ω,
                    dopants :: Vector{GrapheneCoord},
                    τ :: Float64,
                    ϵ :: Float64)
    ΔA = -imag(Δρ_Dopants_Integrand(loc, ω + 1im * η, dopants, τ, ϵ))
    A0 = -imag(Σ(ω + 1im * η))
    return (A0 + ΔA)
end
