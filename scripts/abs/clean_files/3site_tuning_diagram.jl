using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra
using JLD2
using DataFrames
using LaTeXStrings
using CairoMakie
using ChunkSplitters
using ProgressMeter
synceddir(args...) = joinpath(ENV["Dropbox"], "data", "LongerPoorMans", args...)

##
data = []
bdg = true
res = (500, 500)
N = 3
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0, Ez=3)
c = bdg ? FermionBdGBasis(1:N, (:↑, :↓)) : FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
target = LD_cells
f, f!, cache = hamfunc(Hδϕ_Hε(), c, fixedparams)
# εs = sqrt(fixedparams.Ez^2 - first(fixedparams.Δ)^2) .+ range(-0.2, 0.25, length=res[1])
εs = sqrt(fixedparams.Ez^2 - first(fixedparams.Δ)^2) .+ range(-0.15, 0.25, length=res[1]) # 3 sites
δϕs = range(0, pi, length=res[2])
iter = Iterators.product(δϕs, εs) |> collect
hf = δϕε -> f!(cache, δϕε)
eigfunc = δϕε -> diagonalize(f!(cache, δϕε), c)
mapfunc(δϕε) = all_info(eigfunc(δϕε))
##
data = @showprogress map(mapfunc, iter)
##
## Get sweet spots
exps = range(-1.0, 4.0, 6)
prob = ScheduledOptProb(eigfunc, target, GapPenalty(exps))
prob_nodeg = ScheduledOptProb(eigfunc, target)
kwargs = (; iterations=length(exps), initials=[0.5pi, 2.9], ranges=[(0.0, 1.0pi), (2.5, 3.2)])
ss_deg = solve(prob, BestOf(best_algs()); kwargs..., MaxTime=10)
ss_nodeg = solve(prob_nodeg, BestOf(best_algs()); kwargs..., MaxTime=5)
## Save data
wsave(datadir("final_data", "$N-site-tuning.jld2"), Dict("data" => data, "ss_deg" => ss_deg, "ss_nodeg" => ss_nodeg, "εs" => εs, "δϕs" => δϕs, "fixedparams" => fixedparams, "N" => N, "res" => res, "target" => target, "bdg" => bdg))
## Load data
data_dict = load(datadir("final_data", "3-site-tuning.jld2"));
@unpack ss_deg, ss_nodeg, data, εs, δϕs = data_dict;
##
sweet_spots = [[1.4, 2.81], [1.825, 2.9], [2.2, 2.945], ss_deg.sol, ss_nodeg.sol]
sweet_spot_fs = mapfunc.(sweet_spots)
cbwidth = 10
levelcolors = [:darkorange1, :crimson]
linewidth = 1.3
contourcolormap = cgrad(:managua, rev=true)#:vanimo#:managua
contourstyle = :dash
fig = with_theme(theme_latexfonts()) do
    fig = Figure(size=0.75 .* (600, 400), fontsize=20)
    ax = Axis(fig[1, 1]; xlabel=L"ε/Δ", ylabel=L"δϕ", yticks=(pi * [0, 1 / 2, 1], [L"0", L"\frac{\pi}{2}", L"π"]))

    _contourcolormap = cgrad(contourcolormap, rev=true, scale=x -> atanh(1.95 * (x - 1 / 2)))
    colorrange = 0.1 .* (-1, 1)
    levels = 0.1 * range(-1, 1, 15)

    f_econtour = contourf!(ax, εs, δϕs, map(x -> x.gap, data)'; linewidth, levels, linestyle=contourstyle, colormap=_contourcolormap, colorrange)
    hmap = heatmap!(ax, εs, δϕs, map(LD_cells, data)'; colormap=Reverse(:viridis), colorscale=identity, colorrange=(0, 1))

    _f_econtour = contour!(ax, εs, δϕs, map(x -> x.gap, data)'; linewidth, levels, linestyle=contourstyle, colormap=_contourcolormap, colorrange)

    contour!(ax, εs, δϕs, map(x -> x.gap, data)'; linewidth, levels=[0], color=cgrad(contourcolormap, 3, categorical=true)[2])

    markers = [:vline, :cross, :hline, :x][1:3]
    colors = fill(:crimson, 3)
    markersizes = [30, 20, 30]
    f_ss = [scatter!(ax, last(ss), first(ss); marker, markersize, strokewidth=1, color) for (ss, marker, color, markersize) in zip(sweet_spots, markers, colors, markersizes)]

    axislegend(ax, f_ss, ["Phase protection", "\"Topological\"", "Level protection", "Sweet spot"][1:length(f_ss)], position=:rb, labelsize=13, titlesize=13)

    Colorbar(fig[1, 2], hmap; width=cbwidth, ticksize=cbwidth, tickalign=true, ticks=([0, 1 / 2, 1], ["0", "0.5", "1"]), ticklabelsize=15)
    Colorbar(fig[1, 2], hmap; width=cbwidth, ticksize=cbwidth, tickalign=true, ticks=(map(target, sweet_spot_fs[1:3]), ["▌  ", "➕   ", "▬   "]), ticklabelalign=(:right, :center), ticklabelsize=15, ticklabelcolor=:crimson, tickwidth=1, tickcolor=:crimson, ticklabelfont=:bold) # ["ꞁ▮▌  ", "⁺ ➕ ", "╴▬➖ "]
    Colorbar(fig[1, 3], f_econtour; width=cbwidth, ticksize=cbwidth, tickalign=true, ticks=WilkinsonTicks(3), ticklabelsize=15)

    Label(fig[1, 2, Bottom()], " LD", tellwidth=false, tellheight=false, fontsize=20)
    Label(fig[1, 3, Bottom()], "    δE/Δ", tellwidth=false, tellheight=false, fontsize=20)
    colsize!(fig.layout, 3, 0.01)
    colgap!(fig.layout, 2, 10)
    resize_to_layout!(fig)
    fig
end
##
save(plotsdir("sweet_spot_tuning_3_site.png"), fig; px_per_unit=10)

## Calculate perturbative data
# a = FermionBdGBasis(1:3)#; qn=QuantumDots.parity)
a = FermionBasis(1:3; qn=QuantumDots.parity)
# function perturbative_solution(a, M, kwargs...)
#     H = perturbative_hamiltonian(a, M; kwargs...)
#     fullsolve(H, a)
# end
fixedparams_pert = merge(fixedparams[[:t, :θ, :Ez]], (; Δ=fixedparams.Δ * [1, 1, 1]))
pert_solve(a, M, δϕ, ε) = fullsolve((perturbative_hamiltonian(a, M; ε=fill(ε, 3), δϕ=fill(δϕ, 2), fixedparams_pert...)), a)
sol_bdg = pert_solve(a, 1, 1.2, 2.0)
##
pert_datas = [Matrix{Any}(undef, res...), Matrix{Any}(undef, res...)]
n = Threads.nthreads()
##
pert_datas = []
@time foreach(1:2) do M
    # Threads.@threads for (ichunk, inds) in enumerate(chunks(iter; n=n))
    #     _f = δϕε -> pert_solve(a, M, δϕε...)
    #     pert_datas[M][inds] = map(_f, iter[inds])
    # end
    _f = δϕε -> pert_solve(a, M, δϕε...)
    push!(pert_datas, @showprogress "Calculating perturbative solutions for M=$M" map(_f, iter))
end

## Save data
wsave(datadir("final_data", "3-site-perturbative.jld2"), Dict("data" => pert_datas, "εs" => εs, "δϕs" => δϕs, "fixedparams" => fixedparams_pert))
## Laod data
pert_data_dict = load(datadir("final_data", "3-site-perturbative.jld2"));
pert_datas = pert_data_dict["data"];

## Figure
fig_pert = with_theme(theme_latexfonts()) do
    fig = Figure(size=0.8 .* (600, 300), fontsize=20, figure_padding=2)

    g = fig[1, 1] = GridLayout()
    xticks = WilkinsonTicks(3)
    ax = Axis(g[1, 3]; xlabel=L"ε/Δ", ylabel=L"δϕ", yticks=(pi * [0, 1 / 2, 1], [L"0", L"\frac{\pi}{2}", L"π"]), title=L"H", xticks)
    ax1 = Axis(g[1, 1]; title=L"H_1^\text{eff}", xlabel=L"ε/Δ", ylabel=L"δϕ", yticks=(pi * [0, 1 / 2, 1], [L"0", L"\frac{\pi}{2}", L"π"]), xticks)
    ax2 = Axis(g[1, 2]; title=L"H_2^\text{eff}", xlabel=L"ε/Δ", ylabel=L"δϕ", yticks=(pi * [0, 1 / 2, 1], [L"0", L"\frac{\pi}{2}", L"π"]), xticks)
    linkaxes!(ax, ax1, ax2)
    hideydecorations!(ax)
    hideydecorations!(ax2)

    target = x -> MP(x)#x -> norm(x.reduced.two_cells)
    target_label = " MP"    #" LD"
    target_label = "   1-MP"    #" LD"

    hmap_kwargs = (; colormap=Reverse(:viridis), colorscale=log10, colorrange=(1e-3, 1))
    contour_kwargs = (; color=:red, levels=[0.0], linewidth=1.3)

    hmap = heatmap!(ax, εs, δϕs, map(target, data)'; hmap_kwargs...)
    f_cont = contour!(ax, εs, δϕs, map(x -> x.gap, data)'; contour_kwargs...)
    hmap1 = heatmap!(ax1, εs, δϕs, map(target, pert_datas[1])'; hmap_kwargs...)
    f_cont1 = contour!(ax1, εs, δϕs, map(x -> x.gap, pert_datas[1])'; contour_kwargs...)
    hmap2 = heatmap!(ax2, εs, δϕs, map(target, pert_datas[2])'; hmap_kwargs...)
    f_cont2 = contour!(ax2, εs, δϕs, map(x -> x.gap, pert_datas[2])'; contour_kwargs...)

    markers = [:vline, :cross, :hline, :x][1:3]
    colors = fill(:crimson, 3)
    markersizes = [30, 20, 30]

    ticks = ([0, 1 / 2, 1], ["0", "0.5", "1"])
    # ticks=([1e-3,1e-2,1e-1], ["0", "0.5", "1"])
    ticks = LogTicks(WilkinsonTicks(3))
    Colorbar(fig[1, 2], hmap; width=cbwidth, ticksize=cbwidth, tickalign=true, ticks, ticklabelsize=20)

    # Label(fig[1, 2, Top()], target_label; valign=:top, tellheight=false, fontsize=20)
    Label(fig[1, 2, Bottom()], target_label; valign=:center, tellheight=false, fontsize=20)
    # Label(fig[1, 3, Bottom()], "    δE/Δ", tellwidth=false, tellheight=false, fontsize=20)
    # colsize!(fig.layout, 3, 0.01)
    # colgap!(fig.layout, 2, 10)
    # resize_to_layout!(fig)
    fig
end
##
fig_pert_vert = with_theme(theme_latexfonts()) do
    fig = Figure(size=0.75 .* (400, 600), fontsize=20)

    g = fig[1, 1] = GridLayout()
    xticks = WilkinsonTicks(3)

    ax = Axis(g[3, 1]; xlabel=L"ε/Δ", ylabel=L"δϕ", yticks=(pi * [0, 1 / 2, 1], [L"0", L"\frac{\pi}{2}", L"π"]), xticks)
    Label(g[3, 1, TopLeft()], L"H"; halign=:left, tellheight=false)

    ax1 = Axis(g[1, 1]; xlabel=L"ε/Δ", ylabel=L"δϕ", yticks=(pi * [0, 1 / 2, 1], [L"0", L"\frac{\pi}{2}", L"π"]), xticks)
    Label(g[1, 1, TopLeft()], L"H_1^\text{eff}"; halign=:left, tellheight=false)

    ax2 = Axis(g[2, 1]; xlabel=L"ε/Δ", ylabel=L"δϕ", yticks=(pi * [0, 1 / 2, 1], [L"0", L"\frac{\pi}{2}", L"π"]), xticks)
    Label(g[2, 1, TopLeft()], L"H_2^\text{eff}"; halign=:left, tellheight=false)
    linkaxes!(ax, ax1, ax2)
    hidexdecorations!(ax1)
    hidexdecorations!(ax2)

    target = x -> MP(x)#x -> norm(x.reduced.two_cells)
    target_label = " MP"    #" LD"
    target_label = "   1-MP"    #" LD"

    hmap_kwargs = (; colormap=Reverse(:viridis), colorscale=log10, colorrange=(1e-3, 1))
    contour_kwargs = (; color=:red, levels=[0.0], linewidth=1.3)

    hmap = heatmap!(ax, εs, δϕs, map(target, data)'; hmap_kwargs...)
    f_cont = contour!(ax, εs, δϕs, map(x -> x.gap, data)'; contour_kwargs...)
    hmap1 = heatmap!(ax1, εs, δϕs, map(target, pert_datas[1])'; hmap_kwargs...)
    f_cont1 = contour!(ax1, εs, δϕs, map(x -> x.gap, pert_datas[1])'; contour_kwargs...)
    hmap2 = heatmap!(ax2, εs, δϕs, map(target, pert_datas[2])'; hmap_kwargs...)
    f_cont2 = contour!(ax2, εs, δϕs, map(x -> x.gap, pert_datas[2])'; contour_kwargs...)

    markers = [:vline, :cross, :hline, :x][1:3]
    colors = fill(:crimson, 3)
    markersizes = [30, 20, 30]

    ticks = ([0, 1 / 2, 1], ["0", "0.5", "1"])
    # ticks=([1e-3,1e-2,1e-1], ["0", "0.5", "1"])
    ticks = LogTicks(WilkinsonTicks(3))
    Colorbar(fig[1, 2], hmap; width=cbwidth, ticksize=cbwidth, tickalign=true, ticks, ticklabelsize=15)

    # Label(fig[1, 2, Top()], target_label; valign=:top, tellheight=false, fontsize=20)
    Label(fig[1, 2, Bottom()], target_label; valign=:center, tellheight=false, fontsize=20)
    # Label(fig[1, 3, Bottom()], "    δE/Δ", tellwidth=false, tellheight=false, fontsize=20)
    # colsize!(fig.layout, 3, 0.01)
    # colgap!(fig.layout, 2, 10)
    # resize_to_layout!(fig)
    fig
end
##
save(plotsdir("3_site_tuning_perturbative2.png"), fig_pert; px_per_unit=10)
