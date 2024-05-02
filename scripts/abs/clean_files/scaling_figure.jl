using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra
using JLD2
using DataFrames
using LaTeXStrings
using CairoMakie
synceddir(args...) = joinpath(ENV["Dropbox"], "data", "LongerPoorMans", args...)
includet(scriptsdir("abs", "phase_plots.jl"))
includet(scriptsdir("abs", "phase_misc.jl"))
load_full_data(N, folder, dir=synceddir) = wload(dir("phase_diagram", folder, "full_N=$(N)_fixedparams=(t = 0.5, θ = QuantumDots.DiffChainParameter{Float64}(2.746801533890032), V = 0, Δ = 1, U = 0.0, Ez = 3).jld2"))
##
Ns = 2:20
data_ss_deg = Dict(zip(Ns, load_full_data.(Ns, "ss_bdg_noexcgap_bestof2", datadir)));
# data_ss_possible_egap = Dict(zip(Ns, load_full_data.(Ns, "ss_bdg_noexcgap_bestof", synceddir)));
# data_ss_nodeg = Dict(zip(Ns, load_full_data.(Ns, "ss_bdg_noexcgap_bestof_possible_egap", datadir)));
data_ss_nodeg = Dict(zip(Ns, load_full_data.(Ns, "ss_bdg_noexcgap_bestof_nodeg", datadir)));

##
# Ns = sort(filter(k -> k ∉ [], collect(keys(data_ss))))[1:20]
cbwidth = 10
linewidth = 2
markers = [:rect, :circle]
xticks = LinearTicks(5)
xlabel = "Chain length N"
markersize = [14, 12]
LDs_deg = map(N -> LDbdg(data_ss[N]["ss"].optsol), Ns)
LDs_nodeg = map(N -> LDbdg(data_ss_nodeg[N]["ss"].optsol), Ns)
gap_deg = map(N -> abs(data_ss[N]["ss"].optsol.gap), Ns)
gap_nodeg = map(N -> abs(data_ss_nodeg[N]["ss"].optsol.gap), Ns)
excgap_deg = map(N -> data_ss[N]["ss"].optsol.excgap, Ns)
excgap_nodeg = map(N -> data_ss_nodeg[N]["ss"].optsol.excgap, Ns)

fig = with_theme(theme_latexfonts()) do
    fig = Figure(size=0.7 .* (600, 600), fontsize=20)
    ax = Axis(fig[1, 1]; xlabel, ylabel=L"\text{LD-sp}", yscale=log10, yticks=LogTicks(LinearTicks(3)), xticks)
    ylims!(ax, (1e-5, 1))

    f_ld_deg = scatterlines!(ax, Ns, LDs_deg; linewidth, legend=false, markers=true, label="δE ≈ 0", marker=markers[1], markersize=markersize[1])
    f_ld_nodeg = scatterlines!(ax, Ns, LDs_nodeg; linewidth, legend=false, markers=true, label="δE unconstrained", marker=markers[2], markersize=markersize[2])

    axislegend(ax, position=:lb, labelsize=13, titlesize=13)

    ax2 = Axis(fig[2, 1]; xlabel, ylabel="|δE|/Δ", yscale=log10, yticks=LogTicks(LinearTicks(4)), xticks)
    hidexdecorations!(ax; ticks=true, grid=false)
    ylims!(ax2, (1e-18, 1))
    f_gap_deg = scatterlines!(ax2, Ns, gap_deg; linewidth, legend=false, markers=true, label="δE = 0", marker=markers[1], markersize=markersize[1])
    f_gap_nodeg = scatterlines!(ax2, Ns, gap_nodeg; linewidth, legend=false, markers=true, label="δE != 0", marker=markers[2], markersize=markersize[2])

    ax3 = Axis(fig[3, 1]; xlabel, ylabel=L"\text{Excgap}/Δ", yticks=WilkinsonTicks(3), xticks)
    ylims!(ax3, (0, 0.2))
    hidexdecorations!(ax2; ticks=true, grid=false)
    f_excgap_deg = scatterlines!(ax3, Ns, excgap_deg; linewidth, legend=false, markers=true, label="δE = 0", marker=markers[1], markersize=markersize[1])
    f_excgap_nodeg = scatterlines!(ax3, Ns, excgap_nodeg; linewidth, legend=false, markers=true, label="δE != 0", marker=markers[2], markersize=markersize[2])


    # Label(fig[1, 2, Bottom()], " LD", tellwidth=false, tellheight=false, fontsize=20)
    # Label(fig[1, 3, Bottom()], "    δE/Δ", tellwidth=false, tellheight=false, fontsize=20)
    # colsize!(fig.layout, 3, 0.01)
    # colgap!(fig.layout, 2, 10)
    rowgap!(fig.layout, 1, 4)
    rowgap!(fig.layout, 2, 4)
    # resize_to_layout!(fig)
    fig
end
##
save(plotsdir("scaling.png"), fig; px_per_unit=10)