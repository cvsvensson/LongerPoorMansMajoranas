using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra
using JLD2
using DataFrames
using LaTeXStrings
using CairoMakie
## Calculate data
function get_sweet_spots(config)
    @unpack N, fixedparams, exps, optparams, initials, ranges, MaxTime0 = config
    bdg = iszero(fixedparams.U) && iszero(fixedparams.V)
    iterations = length(exps)
    c = bdg ? FermionBdGBasis(1:N, (:↑, :↓)) : FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
    alg = BestOf(best_algs())
    target = LD_cells#bdg ? LDbdg : LD
    f, f!, cache = hamfunc(optparams, c, fixedparams)
    eigfunc = x -> diagonalize(f!(cache, x), c)
    prob = ScheduledOptProb(eigfunc, target, GapPenalty(exps))
    prob_nodeg = ScheduledOptProb(eigfunc, target)
    kwargs = (; iterations, MaxTime=length(alg.optimizers) * MaxTime0 * N, initials, ranges)
    ss = solve(prob, alg; kwargs...)
    ss_nodeg = solve(prob_nodeg, alg; kwargs...)
    @strdict N fixedparams ss ss_nodeg optparams
end
##
ranges = [[(0.0, 1.0pi), (2.5, 3.2)]]
initials = [[1.9, 2.9]]
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0, Ez=3)
exps = range(-1.0, 4.0, 6)
MaxTime0 = 10
optparams = Hδϕ_Hε()
N = collect(2:20)
config = @dict N fixedparams exps optparams initials ranges MaxTime0
configs = dict_list(config)
##
folder = datadir("final_data", "sweet_spot_scaling_ld3")
datas = [produce_or_load(get_sweet_spots, config, folder; filename=x -> savename(x; ignores=("MaxTime0",), allowedtypes=(Int, NamedTuple)))[1] for config in configs];
##
_datas = datas[1:end]
sweet_spots = map(d -> first(filter(x -> abs(x.optsol.gap) < 1e-10, d["ss"].all_sols)), _datas)
sweet_spots_nodeg = map(d -> d["ss_nodeg"], _datas)
##
Ns = map(d -> d["N"], _datas)
LDs_deg = map(d -> LD_cells(d.optsol), sweet_spots)
LDs_nodeg = map(d -> LD_cells(d.optsol), sweet_spots_nodeg)
LDbdgs_deg = map(d -> LDbdg(d.optsol), sweet_spots)
LDbdgs_nodeg = map(d -> LDbdg(d.optsol), sweet_spots_nodeg)
gap_deg = map(d -> abs(d.optsol.gap), sweet_spots)
gap_nodeg = map(d -> abs(d.optsol.gap), sweet_spots_nodeg)
excgap_deg = map(d -> d.optsol.excgap, sweet_spots)
excgap_nodeg = map(d -> d.optsol.excgap, sweet_spots_nodeg)
##
# Ns = sort(filter(k -> k ∉ [], collect(keys(data_ss))))[1:20]
cbwidth = 10
linewidth = 2
markers = [:rect, :circle]
xticks = LinearTicks(5)
xlabel = "Chain length N"
markersize = [14, 12]

fig = with_theme(theme_latexfonts()) do
    fig = Figure(size=0.6 .* (600, 600), fontsize=20)
    ax = Axis(fig[1, 1]; xlabel, ylabel="LD", yscale=log10, yticks=LogTicks(LinearTicks(3)), xticks)
    ylims!(ax, (1e-5, 1))
    c1 = Cycled(1)
    c2 = Cycled(3)
    strokewidth = 1
    common_kwargs = (; linewidth, strokewidth)
    kwargs1 = (; color=c1, marker=markers[1], markersize=markersize[1], label=L"\mathrm{mLD}_0")
    kwargs2 = (; color=c2, marker=markers[2], markersize=markersize[2], label=L"\mathrm{mLD}")
    f_ld_deg = scatterlines!(ax, Ns, LDs_deg; common_kwargs..., kwargs1...)
    f_ld_nodeg = scatterlines!(ax, Ns, LDs_nodeg; common_kwargs..., kwargs2...)
    axislegend(ax, position=:lb, labelsize=13, titlesize=13)


    ax2 = Axis(fig[2, 1]; xlabel, ylabel="|δE|/Δ", yscale=log10, yticks=LogTicks(LinearTicks(4)), xticks)
    hidexdecorations!(ax; ticks=true, grid=false)
    ylims!(ax2, (1e-18, 1))
    f_gap_deg = scatterlines!(ax2, Ns, gap_deg; common_kwargs..., kwargs1...)
    f_gap_nodeg = scatterlines!(ax2, Ns, gap_nodeg; common_kwargs..., kwargs2...)

    ax3 = Axis(fig[3, 1]; xlabel, ylabel=L"E_\text{ex}/Δ", yticks=WilkinsonTicks(3), xticks)
    ylims!(ax3, (0, 0.2))
    hidexdecorations!(ax2; ticks=true, grid=false)
    f_excgap_deg = scatterlines!(ax3, Ns, excgap_deg; common_kwargs..., kwargs1...)
    f_excgap_nodeg = scatterlines!(ax3, Ns, excgap_nodeg; common_kwargs..., kwargs2...)


    Label(fig[1, 1, TopLeft()], "(a)", tellwidth=false, tellheight=false, fontsize=20, padding=(0, 80, 0, 0))
    Label(fig[2, 1, TopLeft()], "(b)", tellwidth=false, tellheight=false, fontsize=20, padding=(0, 80, 0, 0))
    Label(fig[3, 1, TopLeft()], "(c)", tellwidth=false, tellheight=false, fontsize=20, padding=(0, 80, 0, 0))

    rowgap!(fig.layout, 1, 4)
    rowgap!(fig.layout, 2, 4)
    # resize_to_layout!(fig)
    fig
end
##
save(plotsdir("scaling.png"), fig; px_per_unit=10)

##
fig = with_theme(theme_latexfonts()) do
    fig = Figure(size=0.7 .* (600, 400), fontsize=20)
    ax = Axis(fig[1, 1]; xlabel, ylabel="LD", yscale=log10, yticks=LogTicks(LinearTicks(3)), xticks)
    ylims!(ax, (1e-5, 1))
    c1 = Cycled(1)
    c2 = Cycled(3)
    strokewidth = 1
    common_kwargs = (; linewidth, strokewidth)
    kwargs1 = (; color=c1, marker=markers[1], markersize=markersize[1], label="δE ≈ 0")
    kwargs2 = (; color=c2, marker=markers[2], markersize=markersize[2], label="δE unconstrained")
    f_ld_deg = scatterlines!(ax, Ns, LDs_deg; common_kwargs..., kwargs1...)
    fig
end