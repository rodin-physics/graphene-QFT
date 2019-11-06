include("Lib/general.jl")

A_Lattice = readdlm("Data/Density/rho_mu_04_Dop2_A.dat")
B_Lattice = readdlm("Data/Density/rho_mu_04_Dop2_B.dat")

data = Data_Process(A_Lattice, B_Lattice)

XS = data[:, 1]
YS = data[:, 2]
ρ = data[:, 3]

ρ = ρ * 1e4;

bound = 1e-4;
bound = bound * 1e4;

pyplot();
scatter(XS,YS,
        marker_z = ρ,
        markerstrokecolor = :white,
        markerstrokewidth = 0.001,
        markersize = 2.5,
        leg = false,
        aspect_ratio = 1,
        xtickfont = font(12, "Serif"),
        ytickfont = font(12, "Serif"),
        ylims = (-60,60),
        xlims = (-60,60),
        size = (500,400),
        color = :coolwarm,
        clim = (-bound, bound),
        colorbar = true,
        colorbar_title = L"\Delta\rho\times 10^4"
        )
println("Plot Done")
savefig("Test.pdf")
println("Plot Saved")

data
