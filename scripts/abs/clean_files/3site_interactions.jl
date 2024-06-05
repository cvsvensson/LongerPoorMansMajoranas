using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra
using JLD2
using DataFrames
using LaTeXStrings
using CairoMakie
using Accessors
using ProgressMeter
synceddir(args...) = joinpath(ENV["Dropbox"], "data", "LongerPoorMans", args...)

##
res = (5, 5)
N = 3
c = FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
Vs = range(0.0, 0.5, length=res[1])
Us = range(0.0, 8, length=res[2])
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), Δ=1, Ez=3)
target = LD_cells
function get_sweet_spot(U, V, fixedparams=fixedparams)
    fixedparams = @insert fixedparams.U = U
    fixedparams = @insert fixedparams.V = V
    f, f!, cache = hamfunc(Hδϕ_Hε(), c, fixedparams)

    eigfunc = δϕε -> diagonalize(f!(cache, δϕε), c)
    prob = ScheduledOptProb(eigfunc, LD_cells)
    kwargs = (; iterations=4, initials=[0.5pi, 2.9], ranges=[(0.0, 1.0pi), (2 - U / 2 - V, 4 + U / 2 + V)], MaxTime=1, local_reltol=1e-6, verbosity=0, abstol=1e-6)
    exps = range(0.0, 4.0, kwargs.iterations)
    prob_deg = ScheduledOptProb(eigfunc, target, GapPenalty(exps))
    alg = BestOf(reduce(vcat, best_algs()[1:2] for k in 1:1))
    sol_nodeg = solve(prob, alg; kwargs...)
    sol_deg = solve(prob_deg, alg; kwargs...)

    # kwargs_NL = (; length=length(kwargs.initials), tol=1e-10, ftol_rel=kwargs.local_reltol, ftol_abs=kwargs.abstol, xtol_rel=kwargs.local_reltol, xtol_abs=kwargs.abstol, alg = :AUGLAG)
    # prob_NL_deg = LongerPoorMansMajoranas.NLOptProb(eigfunc, target, (sol, x) -> get_gap(sol), kwargs_NL)
    # @time ss_NL_deg = solve(prob_NL_deg; MaxTime=kwargs.MaxTime, initials=kwargs.initials)

    return sol_deg, sol_nodeg
end
##
UViter = Iterators.product(Us, Vs)
data = @showprogress map(UV -> get_sweet_spot(UV...), UViter)
##
## Save data
wsave(datadir("final_data", "UV-tuning_NL.jld2"), Dict("data" => data, "Us" => Us, "Vs" => Vs, "fixedparams" => fixedparams, "N" => N, "res" => res))
## Load data
data_dict = load(datadir("final_data", "3-site-tuning.jld2"));
@unpack data, Us, Vs, fixedparams, N, res = data_dict;
##
cbwidth = 10
levelcolors = [:darkorange1, :crimson]
linewidth = 1.3
contourcolormap = cgrad(:managua, rev=true)#:vanimo#:managua
contourstyle = :dash
cbwidth = 10
linewidth = 1.3
fig_UV = with_theme(theme_latexfonts()) do
    fig = Figure(size=0.8 .* (600, 300), fontsize=20, figure_padding=5)

    g = fig[1, 1] = GridLayout()
    xticks = WilkinsonTicks(3)
    yticks = WilkinsonTicks(3)#(pi * [0, 1 / 2, 1], [L"0", L"\frac{\pi}{2}", L"π"])
    ax = Axis(g[1, 1]; xlabel=L"U_l", title="deg", ylabel=L"U_{nl}", yticks, xticks)
    ax1 = Axis(g[1, 2]; xlabel=L"U_l", title="nodeg", ylabel=L"U_{nl}", yticks, xticks)
    # ax2 = Axis(g[1, 3]; title="δE/Δ", xlabel=L"ε", ylabel=L"δϕ", yticks=(pi * [0, 1 / 2, 1], [L"0", L"\frac{\pi}{2}", L"π"]), xticks)
    linkaxes!(ax, ax1)#, ax2)
    # linkaxes!(ax, ax1, ax2)
    hideydecorations!(ax1)
    # hideydecorations!(ax2)

    target = x -> LD_cells(x)#x -> norm(x.reduced.two_cells)

    contour_kwargs = (; color=:red, levels=[0.0], linewidth=1.3)
    LD_deg = reshape(map(ss -> LD_cells(ss[1].optsol), data), res...)
    LD_nodeg = reshape(map(ss -> LD_cells(ss[2].optsol), data), res...)
    hmap_kwargs = (; colormap=Reverse(:viridis), colorscale=identity, colorrange=(0, maximum([LD_deg..., LD_nodeg...])))
    hmap = heatmap!(ax, Us, Vs, LD_deg; hmap_kwargs...)
    f_egap = heatmap!(ax1, Us, Vs, LD_nodeg; hmap_kwargs...)
    # f_Q = contour!(ax1, εs_high_res, δϕs_high_res, dataQ2; contour_kwargs...)
    # l_Q = lines!(ax1, Float64[], Float64[]; label="Q=0", contour_kwargs...)
    # axislegend(ax1; position=:lt)

    ticks = ([0, 0.25, 1 / 2, 1], ["0", ".25", "0.5", "1"])
    # ticks=([1e-3,1e-2,1e-1], ["0", "0.5", "1"])
    ticklabelsize = 16
    Colorbar(g[1, 3], hmap; width=cbwidth, ticksize=cbwidth, tickalign=true, ticks=LinearTicks(3), ticklabelsize)
    # Label(g[1, 2, Bottom()], " LD", tellwidth=false, tellheight=false, fontsize=20)
    Label(g[1, 3, Bottom()], " LD", tellwidth=false, tellheight=false, fontsize=20)

    # Colorbar(g[1, 3], f_egap; width=cbwidth, ticksize=cbwidth, tickalign=true, ticks=LogTicks(WilkinsonTicks(3)), ticklabelsize)
    # colgap!(g, 1, 10)
    # colgap!(g, 2, 10)
    fig
end

##
heatmap(Us, Vs, reshape(map(ss -> log(abs((ss[1].optsol.gap)/ (ss[2].optsol.gap))), data), res...), colorrange=3 .* (-1, 1), colormap=:redsblues)
##
save(plotsdir("3_site_UV_sweet_spot.png"), fig_UV; px_per_unit=10)
