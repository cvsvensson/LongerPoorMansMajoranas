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
res = (400, 400)
N = 2
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0, Ez=3)
bdg = iszero(fixedparams.U) && iszero(fixedparams.V)
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
prob_level = ScheduledOptProb(eigfunc, target, GapPenalty(exps) + ScheduledPenalty((sol, x, i) -> 10^(-i)abs(get_gap(eigfunc(x .+ [0, 0.01])) - get_gap(sol)) / 0.01))
prob_phase = ScheduledOptProb(eigfunc, target, GapPenalty(exps) + ScheduledPenalty((sol, x, i) -> 10^(-i)abs(get_gap(eigfunc(x .+ [0.01, 0])) - get_gap(sol)) / 0.01))
prob_nodeg = ScheduledOptProb(eigfunc, target)
kwargs = (; iterations=length(exps), initials=[0.5pi, 2.9], ranges=[(0.0, 1.0pi), (2.5, 3.2)], MaxTime=5)
ss_deg = solve(prob, BestOf(best_algs()); kwargs...)
ss_phase = solve(prob, BestOf(best_algs()); kwargs..., ranges=[(0.0, 3pi / 4), (2.7, 2.9)])
ss_level = solve(prob, BestOf(best_algs()); kwargs..., ranges=[(0.0, 1.0pi), (2.9, 3.2)])
ss_nodeg = solve(prob_nodeg, BestOf(best_algs()); kwargs...)
## Save data
wsave(datadir("final_data", "$N-site-tuning.jld2"), Dict("data" => data, "ss_deg" => ss_deg, "ss_nodeg" => ss_nodeg, "εs" => εs, "δϕs" => δϕs, "fixedparams" => fixedparams, "N" => N, "res" => res, "target" => target, "bdg" => bdg))
## Load data
data_dict = load(datadir("final_data", "$N-site-tuning.jld2"));
@unpack ss_deg, ss_nodeg, data, εs, δϕs = data_dict;
##
# sweet_spots = [[1.4, 2.81], [1.825, 2.9], [2.2, 2.945], ss_deg.sol, ss_nodeg.sol]
sweet_spots = [ss_phase.sol, ss_nodeg.sol, ss_level.sol]
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
    levels = 0.2 * range(-1, 1, 15)

    f_econtour = contourf!(ax, εs, δϕs, map(x -> x.gap, data)'; linewidth, levels, linestyle=contourstyle, colormap=_contourcolormap, colorrange)
    hmap = heatmap!(ax, εs, δϕs, map(LD_cells, data)'; colormap=Reverse(:viridis), colorscale=identity, colorrange=(0, 1))

    _f_econtour = contour!(ax, εs, δϕs, map(x -> x.gap, data)'; linewidth, levels, linestyle=contourstyle, colormap=_contourcolormap, colorrange)

    contour!(ax, εs, δϕs, map(x -> x.gap, data)'; linewidth, levels=[0], color=cgrad(contourcolormap, 3, categorical=true)[2])

    markers = [:circle]
    colors = fill(:crimson, 3)
    markersizes = [20, 20, 30]
    f_ss = [scatter!(ax, last(ss), first(ss); marker, markersize, strokewidth=2, color) for (ss, marker, color, markersize) in zip(sweet_spots[1:1], markers, colors, markersizes)]

    #axislegend(ax, f_ss, ["Phase protection", "\"Topological\"", "Level protection", "Sweet spot"][1:length(f_ss)], position=:rb, labelsize=13, titlesize=13)

    Colorbar(fig[1, 2], hmap; width=cbwidth, ticksize=cbwidth, tickalign=true, ticks=([0, 1 / 2, 1], ["0", "0.5", "1"]), ticklabelsize=15)
    #Colorbar(fig[1, 2], hmap; width=cbwidth, ticksize=cbwidth, tickalign=true, ticks=(map(target, sweet_spot_fs[1:3]), ["▌  ", "➕   ", "▬   "]), ticklabelalign=(:right, :center), ticklabelsize=15, ticklabelcolor=:crimson, tickwidth=1, tickcolor=:crimson, ticklabelfont=:bold) # ["ꞁ▮▌  ", "⁺ ➕ ", "╴▬➖ "]
    Colorbar(fig[1, 3], f_econtour; width=cbwidth, ticksize=cbwidth, tickalign=true, ticks=WilkinsonTicks(3), ticklabelsize=15)

    Label(fig[1, 2, Bottom()], " LD", tellwidth=false, tellheight=false, fontsize=20)
    Label(fig[1, 3, Bottom()], "    δE/Δ", tellwidth=false, tellheight=false, fontsize=20)
    colsize!(fig.layout, 3, 0.01)
    colgap!(fig.layout, 2, 10)
    resize_to_layout!(fig)
    fig
end
##
save(plotsdir("sweet_spot_tuning_2_site.png"), fig; px_per_unit=10)
