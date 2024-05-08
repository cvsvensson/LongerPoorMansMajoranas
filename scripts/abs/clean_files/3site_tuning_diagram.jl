using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra#, BlackBoxOptim
using JLD2
using DataFrames
using LaTeXStrings
using CairoMakie
using ChunkSplitters
synceddir(args...) = joinpath(ENV["Dropbox"], "data", "LongerPoorMans", args...)

##
data = []
bdg = false
res = (500, 500)
N = 3
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0, Ez=3)
c = bdg ? FermionBdGBasis(1:N, (:↑, :↓)) : FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
target = bdg ? LDbdg : LD
f, f!, cache = hamfunc(Hδϕ_Hε(), c, fixedparams)
εs = sqrt(fixedparams.Ez^2 - fixedparams.Δ^2) .+ range(-0.15, 0.25, length=res[1])
δϕs = range(0, pi, length=res[2])
iter = Iterators.product(δϕs, εs) |> collect
hf = δϕε -> f!(cache, δϕε)
mapfunc(δϕε) = fullsolve(hf(δϕε), c)
# mapgrad(δϕε) = get_gap_gradient(hf, c, [δϕε...])
# @time grads = Folds.map(mapgrad, iter);
##
data = Matrix{Any}(undef, res...)
n = Threads.nthreads()
caches = [deepcopy(cache) for _ in 1:n]
@time begin
    Threads.@threads for (ichunk, inds) in enumerate(chunks(iter; n=n))
        _f = δϕε -> fullsolve(f!(caches[ichunk], δϕε), c)
        data[inds] = map(_f, iter[inds])
    end
end
# @time data = map(mapfunc, iter);

exps = range(0, 4, 6)
prob = ScheduledOptProb(mapfunc, target, GapPenalty(exps))
prob_nodeg = ScheduledOptProb(mapfunc, target)
ss_sol = solve(prob, BestOf(best_algs()); MaxTime=10)
ss_sol_nodeg = solve(prob, BestOf(best_algs()); MaxTime=5)
##
wsave(datadir("final_data", "3-site-tuning.jld2"), Dict("data" => data, "ss_deg" => ss_sol, "ss_nodeg" => ss_sol_nodeg))
##
sweet_spots = [[2.81, 1.4], [2.9, 1.825], [2.945, 2.2], reverse(ss_sol.sol), reverse(ss_sol2.sol)]
sweet_spot_fs = mapfunc.(reverse.(sweet_spots))
cbwidth = 10
levelcolors = [:darkorange1, :crimson]
linewidth = 1.3
contourcolormap = cgrad(:managua, rev=true)#:vanimo#:managua
contourstyle = :dash
fig = with_theme(theme_latexfonts()) do
    fig = Figure(size=0.75 .* (600, 400), fontsize=20)
    ax = Axis(fig[1, 1]; xlabel=L"ε", ylabel=L"δϕ", yticks=(pi * [0, 1 / 2, 1], [L"0", L"\frac{\pi}{2}", L"π"]))


    _contourcolormap = cgrad(contourcolormap, rev=true, scale=x -> atanh(1.95 * (x - 1 / 2)))
    colorrange = 0.1 .* (-1, 1)
    levels = 0.1 * range(-1, 1, 15)

    f_econtour = CairoMakie.contourf!(ax, εs, δϕs, map(x -> x.gap, data)'; linewidth, levels, linestyle=contourstyle, colormap=_contourcolormap, colorrange)
    hmap = CairoMakie.heatmap!(ax, εs, δϕs, map(target, data)'; colormap=Reverse(:viridis), colorscale=identity, colorrange=(0, 1))

    _f_econtour = CairoMakie.contour!(ax, εs, δϕs, map(x -> x.gap, data)'; linewidth, levels, linestyle=contourstyle, colormap=_contourcolormap, colorrange)

    CairoMakie.contour!(ax, εs, δϕs, map(x -> x.gap, data)'; linewidth, levels=[0], color=cgrad(contourcolormap, 3, categorical=true)[2])

    # f_E0 = lines!(ax, Float64[], Float64[]; linewidth, color=levelcolors[2])
    # f_Es = lines!(ax, Float64[], Float64[]; linewidth, linestyle=:dash, color=levelcolors[1])

    markers = [:vline, :cross, :hline, :x][1:3]
    colors = fill(:crimson, 3)#[:red, :red, :red]#[:aqua, :red, :aqua]
    markersizes = [30, 20, 30]
    f_ss = [CairoMakie.scatter!(ax, first(ss), last(ss); marker, markersize, strokewidth=1, color) for (ss, marker, color, markersize) in zip(sweet_spots, markers, colors, markersizes)]

    #axislegend(ax, [f_E0, f_Es], ["δE = 0", "δE = 0.01Δ"], position=:lb, labelsize=17)
    axislegend(ax, f_ss, ["Phase protection", "\"Topological\"", "Level protection", "Sweet spot"][1:length(f_ss)], position=:rb, labelsize=13, titlesize=13)
    # Colorbar(fig[1, 2], hmap; width=cbwidth, ticksize=cbwidth, tickalign=true, ticks=[0, 1 // 2, 1], ticklabelsize=15)

    Colorbar(fig[1, 2], hmap; width=cbwidth, ticksize=cbwidth, tickalign=true, ticks=([0, 1 / 2, 1], ["0", "0.5", "1"]), ticklabelsize=15)
    Colorbar(fig[1, 2], hmap; width=cbwidth, ticksize=cbwidth, tickalign=true, ticks=(map(target, sweet_spot_fs[1:3]), ["▌  ", "➕   ", "▬   "]), ticklabelalign=(:right, :center), ticklabelsize=15, ticklabelcolor=:crimson, tickwidth=1, tickcolor=:crimson, ticklabelfont=:bold) # ["ꞁ▮▌  ", "⁺ ➕ ", "╴▬➖ "]
    Colorbar(fig[1, 3], f_econtour; width=cbwidth, ticksize=cbwidth, tickalign=true, ticks=WilkinsonTicks(3), ticklabelsize=15)

    # Label(fig[1, 2, Bottom()], MarkerElement(color=:blue, marker='x', markersize=15,
    # strokecolor=:black), tellwidth=false, tellheight=false, fontsize=20)
    Label(fig[1, 2, Bottom()], " LD", tellwidth=false, tellheight=false, fontsize=20)
    Label(fig[1, 3, Bottom()], "    δE/Δ", tellwidth=false, tellheight=false, fontsize=20)
    colsize!(fig.layout, 3, 0.01)
    colgap!(fig.layout, 2, 10)
    resize_to_layout!(fig)
    fig
end
##
save(plotsdir("sweet_spot_tuning_3_site.png"), fig; px_per_unit=10)
#colsize!(fig.layout, 1, Aspect(1, 1.0));
#colgap!(fig.layout, 7);

##
fig = Figure(size=(900, 400), fontsize=20)
axs = [Axis(fig[1, j], aspect=1) for j in 1:2]
cmap = :roma
contour!(axs[1], x, y, fargs, levels=30, colormap=cmap)
pltobj1 = heatmap!(axs[2], x, y, fargs, colorrange=(-π, π), colormap=cmap)
contour!(axs[2], x, y, fargs, levels=30, color=:white, linewidth=0.85)
Colorbar(fig[1, 3], pltobj1, scale=:log10)
limits!(axs[1], -2, 2, -2, 2)
fig


plot_f(data2, target; clim_max=1, leg=:bottomleft, plot_ss=false, title="3 sites", colorbar_title="LD", ss_label=false)
scatter!(first.(sweet_spots), last.(sweet_spots))