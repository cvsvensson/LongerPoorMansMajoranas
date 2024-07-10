using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra
using CairoMakie
using JLD2
using DataFrames
using LaTeXStrings
using ProgressMeter
##
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0, Ez=3)
bdg = iszero(fixedparams.V) && iszero(fixedparams.U)
target = LD_cells
## 2 site sweet  spot
c2 = bdg ? FermionBdGBasis(1:2, (:↑, :↓)) : FermionBasis(1:2, (:↑, :↓); qn=QuantumDots.parity)
f2, f2!, cache2 = hamfunc(Hδϕ_Hε(), c2, fixedparams)
exps = range(-1, 5, 6)
eigfunc = x -> diagonalize(f2!(cache2, x), c2)
prob_2_site = ScheduledOptProb(eigfunc, target, GapPenalty(exps))
ss2 = solve(prob_2_site, BestOf(best_algs()); iterations=length(exps), MaxTime=1, initials=[2, 2], ranges=[(0.0, 1.0pi), (0.0, 5.0)])

## 3 site sweet spot
initials = [2.2, sqrt(fixedparams.Ez^2 - first(fixedparams.Δ)^2)]
c = bdg ? FermionBdGBasis(1:3, (:↑, :↓)) : FermionBasis(1:3, (:↑, :↓); qn=QuantumDots.parity)
f, f!, cache = hamfunc(Rδϕ_Rε(), c, fixedparams);
fh, fh!, cacheh = hamfunc(Hδϕ_Hε(), c, fixedparams);
@assert norm(f([1, 1, 1]) - fh([1, 1])) < 1e-16

eigfunc = x -> diagonalize(fh!(cacheh, x), c)
prob_3_site_homogeneous = ScheduledOptProb(eigfunc, target, GapPenalty(exps))
ss3_homogeneous = solve(prob_3_site_homogeneous, BestOf(best_algs()); iterations=length(exps), MaxTime=2, initials=[2, 2], ranges=[(0.0, 1.0pi), (0.0, 5.0)])
##
ss = ss3_homogeneous
ε0 = ss.sol[2]
initials = collect(ss.sol)
MaxTime = 10
minexcgap = ss.optsol.excgap * 0.0
phase_data = []
level_data = []
level_data_pf = []
phase_data_pf = []
nodeg_data = []
δε2s = range(-0.2, 0.45, length=15)
@showprogress for δε2 in δε2s
    hamfunc = δϕε1 -> f!(cache, [δϕε1[1], δϕε1[2], δϕε1[2] + δε2])
    eigfunc = x -> diagonalize(hamfunc(x), c)
    prob = ScheduledOptProb(eigfunc, target, GapPenalty(exps))
    prob_nodeg = ScheduledOptProb(eigfunc, target)
    #WARNING: These ranges are a hack to get the different branches for these particular parameters
    phase_ranges = [(0.0, 2.0), 2.82 .+ 1 .* (-2, 0.0)] #phase branch
    level_ranges = [(0.0, 1.0pi), 2.85 .+ 1 .* (0.0, 0.3)] # level branch
    nodeg_ranges = [(0.0, 1.0pi), ε0 .+ 1 .* (-0.3, 0.3)]

    iterations = length(exps)
    kwargs = (; iterations, MaxTime, initials, verbosity=0)
    # phase branch
    phase_sol = solve(prob, BestOf(best_algs()); ranges=phase_ranges, kwargs...)
    push!(phase_data, (phase_sol))
    # level branch
    level_sol = solve(prob, BestOf(best_algs()); ranges=level_ranges, kwargs...)
    push!(level_data, (level_sol))
    # non-degenerate branch
    nodeg_sol = solve(prob_nodeg, BestOf(best_algs()); kwargs..., MaxTime=MaxTime / 3, initials=[2.0, ε0 + 0.1], ranges=nodeg_ranges)
    push!(nodeg_data, (nodeg_sol))
end
##
wsave(datadir("final_data", "branches", "$N-site-branches.jld2"), Dict("nodeg_data" => nodeg_data, "phase_data" => phase_data, "level_data" => level_data, "δε2s" => δε2s, "fixedparams" => fixedparams, "N" => N, "ss2" => ss2, "ss3_homogeneous" => ss3_homogeneous))
## Load data
data_dict = load(datadir("final_data", "branches", "3-site.jld2"))
@unpack nodeg_data, phase_data, level_data, δε2s, fixedparams, N, ss2, ss3_homogeneous = data_dict
##
function branch_plot!(ax, δε2s, datas, two_site, target)
    strokewidth = 1
    markersize = 10
    markers = [:circle, :rect, :(+)]
    # phase_marker = :circle#:vline #'▌'
    # level_marker = :rect#:hline #'▬'
    targets = map(data -> map(target, data), datas)
    colors = Makie.wong_colors()[[1, 3, 4]]

    figs = [CairoMakie.scatterlines!(ax, δε2s, data; color=color, marker=marker, strokewidth, markersize) for (data, color, marker) in zip(targets, colors, markers)]
    f_2s = hlines!(ax, [target(two_site)], label="2 site sweet spot", linewidth=2, linestyle=:dash, color=:black)
    # f_homogeneous = vlines!(ax, [0]; color=:red, linestyle=:dot)

    return figs, f_2s#, f_homogeneous
end
function select_degenerate_sweet_spot(data, cutoff)
    i = findfirst(x -> abs(x.optsol.gap) < cutoff, data.all_sols)
    data.all_sols[i]
end

curve_labels = ["Phase branch", "Level branch", L"\mathrm{mLD} \text{ (3-site)}", L"\mathrm{mLD}_0 \text{ (2-site)}"]#, "Homogeneous"]
curve_labels = [L"\mathrm{mLD}_0 \text{ (Phase branch)}", L"\mathrm{mLD}_0 \text{ (Level branch)}", L"\mathrm{mLD}", L"\mathrm{mLD}_0 \text{ (2-site)}"]#, "Homogeneous"]
branch_fig2 = let level_data = level_data, phase_data = phase_data

    level_data = map(level_data) do d
        select_degenerate_sweet_spot(d, 1e-4)
    end
    phase_data = map(phase_data) do d
        select_degenerate_sweet_spot(d, 1e-4)
    end
    with_theme(theme_latexfonts()) do
        fig = Figure(size=0.75 .* (800, 500), fontsize=20)

        gl = fig[1:2, 1] = GridLayout()
        gr = fig[1:2, 2] = GridLayout()

        xticks = LinearTicks(4)
        ax1 = Axis(gl[1, 1]; xlabel=L"δε_2/Δ", ylabel="||∇δE||")
        ax2 = Axis(gl[2, 1]; xlabel=L"δε_2/Δ", ylabel="δϕ", yticks=(pi * [0, 1 / 2, 1], ["0", L"\frac{\pi}{2}", "π"]), xticks)
        ax3 = Axis(gr[2, 1]; xlabel=L"δε_2/Δ", ylabel=L"ε_{13}", xticks)
        linkxaxes!(ax1, ax2, ax3)
        hidexdecorations!(ax1; ticks=true, grid=false)
        CairoMakie.ylims!(ax1, (0, 0.75))
        CairoMakie.ylims!(ax2, (0, pi))
        datas = [phase_data, level_data, nodeg_data]
        fs_ld, f_ld_2s = branch_plot!(ax1, δε2s, datas, ss2, x -> norm(x.gradient))
        branch_plot!(ax2, δε2s, datas, ss2, x -> (x.sol[1]))
        branch_plot!(ax3, δε2s, datas, ss2, x -> (x.sol[2]))

        Label(gl[1, 1, TopLeft()], "(a)", tellwidth=false, tellheight=false, fontsize=20, padding=(0, 60, 0, 0))
        Label(gl[2, 1, TopLeft()], "(b)", tellwidth=false, tellheight=false, fontsize=20, padding=(0, 60, 0, 0))
        Label(gr[2, 1, TopLeft()], "(c)", tellwidth=false, tellheight=false, fontsize=20, padding=(0, 60, 0, 0))

        leg = Legend(gr[1, 1],
            [fs_ld..., f_ld_2s],
            curve_labels)
        leg.tellwidth = false
        fig
    end
end
branch_fig3 = let level_data = level_data, phase_data = phase_data
    with_theme(theme_latexfonts()) do
        fig = Figure(size=0.75 .* (800, 500), fontsize=20)
        gl = fig[1:2, 1] = GridLayout()
        gr = fig[1:2, 2] = GridLayout()
        xticks = LinearTicks(4)
        level_data = map(level_data) do d
            d.all_sols[findfirst(sol -> abs(sol.optsol.gap) < 1e-4, d.all_sols)]
        end
        phase_data = map(phase_data) do d
            d.all_sols[findfirst(sol -> abs(sol.optsol.gap) < 1e-4, d.all_sols)]
        end

        ax1 = Axis(gl[1, 1]; xlabel=L"δε_2/Δ", ylabel="LD", xticks)
        ax2 = Axis(gl[2, 1]; xlabel=L"δε_2/Δ", ylabel=L"E_\text{ex}/\Delta", xticks)
        hidexdecorations!(ax1; ticks=true, grid=false)
        CairoMakie.ylims!(ax1, (0, 0.75))
        ax3 = Axis(gr[2, 1]; xlabel=L"δε_2/Δ", ylabel="δE/Δ", xticks)
        linkxaxes!(ax1, ax2, ax3)
        datas = [phase_data, level_data, nodeg_data]
        fs, f_2s = branch_plot!(ax1, δε2s, datas, ss2, x -> LD_cells(x.optsol))
        branch_plot!(ax2, δε2s, datas, ss2, x -> (x.optsol.excgap))
        branch_plot!(ax3, δε2s, datas, ss2, x -> (x.optsol.gap))

        Label(gl[1, 1, TopLeft()], "(a)", tellwidth=false, tellheight=false, fontsize=20, padding=(0, 60, 0, 0))
        Label(gl[2, 1, TopLeft()], "(b)", tellwidth=false, tellheight=false, fontsize=20, padding=(0, 60, 0, 0))
        Label(gr[2, 1, TopLeft()], "(c)", tellwidth=false, tellheight=false, fontsize=20, padding=(0, 60, 0, 0))

        leg = Legend(gr[1, 1],
            [fs..., f_2s],
            curve_labels)
        leg.tellwidth = false

        fig
    end
end
##
CairoMakie.save(plotsdir("sweet_spot_branches.png"), branch_fig3; px_per_unit=10)
##
CairoMakie.save(plotsdir("sweet_spot_branches_appendix.png"), branch_fig2; px_per_unit=10)
