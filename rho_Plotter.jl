include("Lib/general.jl")

A_Lattice = readdlm("Data/Density/rho_mu_04_Dop2_A.txt")
B_Lattice = readdlm("Data/Density/rho_mu_04_Dop2_B.txt")

data = Data_Process(A_Lattice, B_Lattice)

XS = [data[1][1, :] data[2][1, :]]
YS = [data[1][2, :] data[2][2, :]]
ρ = [data[1][3, :] data[2][3, :]]

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
        # xaxis = (L"\AA", font(20, "Serif")),
        # yaxis = (L"\AA", font(20, "Serif")),
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
# savefig("rho_mu_04_Dop2.pdf")
println("Plot Saved")
