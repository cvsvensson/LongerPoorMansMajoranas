using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra#, BlackBoxOptim
using Symbolics
using Folds
using Accessors
using JLD2
using DataFrames
using LaTeXStrings
using FiniteDiff
using CairoMakie
synceddir(args...) = joinpath(ENV["Dropbox"], "data", "LongerPoorMans", args...)
includet(scriptsdir("abs", "phase_plots.jl"))
includet(scriptsdir("abs", "phase_misc.jl"))

##
data = []
bdg = true
res = (200, 107)
N = 3
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0.0, Ez=3)
c = bdg ? FermionBdGBasis(1:N, (:↑, :↓)) : FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
f, f!, cache = hamfunc(Hδϕ_Hε(), c, fixedparams)
εs = sqrt(fixedparams.Ez^2 - fixedparams.Δ^2) .+ range(-0.15, 0.25, length=res[1])
δϕs = range(0, pi, length=res[2])
iter = Iterators.product(δϕs, εs) |> collect
hf = f#δϕε -> f!(cache, δϕε)
mapfunc(δϕε) = fullsolve(hf(δϕε), c)
mapgrad(δϕε) = get_gap_gradient(hf, c, [δϕε...])
@time data = Folds.map(mapfunc, iter);
@time grads = Folds.map(mapgrad, iter);
data2 = Dict("data" => data, "x" => εs, "y" => δϕs, "labels" => ("ε", "δϕ"), "N" => N)
##
findmin(LDbdg, data)
findmin(x -> abs(x[1]), grads)
findmin(x -> abs(x[1][1]) + 100abs(x[2].gap), zip(grads, data) |> collect)
findmin(x -> abs(x[2]), grads)
##
Theme(
    Lines=(
        linewidth=4,
        linestyle=:dash,
    )
)
##
heatmap(εs, δϕs, map(LDbdg, data), c=:viridis, clim=(0, 1), xlabel="ε", ylabel="δϕ", title="LD")
heatmap(εs, δϕs, map(x -> x[1], grads), c=:redsblues, clim=0.1 .* (-1, 1), xlabel="ε", ylabel="δϕ", title="phase grad")
heatmap(εs, δϕs, map(x -> x[2], grads), c=:redsblues, clim=0.1 .* (-1, 1), xlabel="ε", ylabel="δϕ", title="eps grad")
heatmap(εs, δϕs, map(norm, grads), c=:viridis, clim=(0, 1), xlabel="ε", ylabel="δϕ", title="|grad|")
##
sweet_spots = [[2.81, 1.4], [2.9, 1.82], [2.945, 2.2]]
fig = Figure(size=0.7 .* (600, 400));
ax = Axis(fig[1, 1]; xlabel="ε", ylabel="δϕ");
hmap = heatmap!(ax, εs, δϕs, map(LDbdg, data)'; colormap=Reverse(:viridis), colorscale=log10);
Colorbar(fig[1, 2], hmap; label="LD", width=15, ticksize=15, tickalign=true);
contour!(ax, εs, δϕs, map(x -> x.gap, data)', levels=[-0.01, 0.01], color=:orange)
contour!(ax, εs, δϕs, map(x -> x.gap, data)', levels=[0], color=:red)
scatter!(ax, first.(sweet_spots), last.(sweet_spots), markersize=10, strokewidth=2, color=:red)
fig

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


plot_f(data2, LDbdg; clim_max=1, leg=:bottomleft, plot_ss=false, title="3 sites", colorbar_title="LD", ss_label=false)
scatter!(first.(sweet_spots), last.(sweet_spots))