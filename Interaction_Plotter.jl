include("Lib/general.jl")

A_Lattice = readdlm("Data/Interaction/F_A_mu_02_A.txt")
B_Lattice = readdlm("Data/Interaction/F_A_mu_02_B.txt")

data = Data_Process(A_Lattice, B_Lattice)

XS = [data[1][1, :] data[2][1, :]]
YS = [data[1][2, :] data[2][2, :]]
FI = [data[1][3, :] data[2][3, :]]

bound = 0.2;

bound = bound * 1e0;
#
FI = FI * 1e0

pyplot();
scatter(XS,YS,
        marker_z = FI,
        markerstrokecolor = :white,
        markerstrokewidth = 0.001,
        markersize = 12,
        leg = false,
        aspect_ratio=1,
        # xaxis = (L"\AA", font(20, "Serif")),
        # yaxis = (L"\AA", font(20, "Serif")),
        xtickfont = font(12, "Serif"),
        ytickfont = font(12, "Serif"),
        ylims = (-10,10),
        xlims = (-10,10),
        size = (500,400),
        color = :coolwarm,
        clim = (-bound, bound),
        colorbar = true,
        colorbar_title = L"F_I\,(\mathrm{eV})"
        )
println("Plot Done")
savefig("Test.pdf")
# savefig("F_I_mu_02.pdf")
println("Plot Saved")
