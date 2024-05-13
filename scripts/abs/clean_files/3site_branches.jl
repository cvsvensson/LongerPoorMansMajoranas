using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra
using CairoMakie
using JLD2
using DataFrames
using LaTeXStrings
synceddir(args...) = joinpath(ENV["Dropbox"], "data", "LongerPoorMans", args...)
##
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0, Ez=3)
bdg = iszero(fixedparams.V) && iszero(fixedparams.U)
target = bdg ? LDbdg : LD
## 2 site sweet  spot
c2 = bdg ? FermionBdGBasis(1:2, (:↑, :↓)) : FermionBasis(1:2, (:↑, :↓); qn=QuantumDots.parity)
f2, f2!, cache2 = hamfunc(Hδϕ_Hε(), c2, fixedparams)
exps = range(0, 4, 5)
fs = x -> fullsolve(f2!(cache2, x), c2)
prob_2_site = ScheduledOptProb(fs, target, GapPenalty(exps))
ss2 = solve(prob_2_site, BestOf(best_algs()); iterations=length(exps), maxiters=100000, MaxTime=1, initials=[2, 2], ranges=[(0.0, 1.0pi), (-5.0, 5.0)])

## 3 site sweet spot
initials = [2.2, sqrt(fixedparams.Ez^2 - fixedparams.Δ^2)]#2.95]
c = bdg ? FermionBdGBasis(1:3, (:↑, :↓)) : FermionBasis(1:3, (:↑, :↓); qn=QuantumDots.parity)
f, f!, cache = hamfunc(Rδϕ_Rε(), c, fixedparams)
fh, fh!, cacheh = hamfunc(Hδϕ_Hε(), c, fixedparams)
@assert norm(f([1, 1, 1]) - fh([1, 1])) < 1e-16

fhs = x -> fullsolve(fh!(cacheh, x), c)
prob_3_site_homogeneous = ScheduledOptProb(fhs, target, GapPenalty(exps))
ss3_homogeneous = solve(prob_2_site, BestOf(best_algs()); iterations=length(exps), maxiters=100000, MaxTime=2, initials=[2, 2], ranges=[(0.0, 1.0pi), (0.0, 5.0)])
##
ss = ss3_homogeneous
ε0 = ss.sol[2]
initials = collect(ss.sol)
MaxTime = 1
minexcgap = ss.optsol.excgap * 0.0
phase_data = []
level_data = []
level_data_pf = []
phase_data_pf = []
nodeg_data = []
δε2s = range(-0.2, 0.45, length=15) #Phase branch
# δε2s = range(-0.2, 0.7, length=10) #level branch
for δε2 in δε2s
    # alg = BestOf(best_algs())
    hamfunc = δϕε1 -> f!(cache, [δϕε1[1], δϕε1[2], δϕε1[2] + δε2])
    fs = x -> fullsolve(hamfunc(x), c)
    prob = ScheduledOptProb(fs, target, GapPenalty(exps))
    prob_nodeg = ScheduledOptProb(fs, target)
    #WARNING: This is a hack to get the phase branch
    phase_ranges = [(0.0, 2.0), 2.82 .+ 1 .* (-2, 0.0)] #phase branch
    level_ranges = [(0.0, 1.0pi), 2.85 .+ 1 .* (0.0, 0.3)] # level branch
    nodeg_ranges = [(0.0, 1.0pi), ε0 .+ 1 .* (-0.3, 0.3)]
    # ranges = [(0.0, 1.0pi), ε0 .+ 2 .* (-1, 1)]

    # level_penalty(pf) = ScheduledPenalty((sol, x, i) -> (i in eachindex(pf) ? pf[i] : last(pf)) * abs(get_gap_gradient(fs, x)[2]))
    # phase_penalty(pf) = ScheduledPenalty((sol, x, i) -> (i in eachindex(pf) ? pf[i] : last(pf)) * abs(get_gap_gradient(fs, x)[1]))
    # lp = level_penalty([10, 1, 0.1, 0.01, 0.0001])
    # ϕp = phase_penalty([10, 1, 0.1, 0.01])
    # prob_level_pf = ScheduledOptProb(fs, target, lp + GapPenalty(exps))
    # prob_phase_pf = ScheduledOptProb(fs, target, ϕp + GapPenalty(exps))

    iterations = length(exps)
    kwargs = (; iterations, MaxTime, initials)
    # phase branch
    phase_sol = solve(prob, BestOf(best_algs()); ranges=phase_ranges, kwargs...)
    push!(phase_data, phase_sol)
    # phase branch with penalty
    # phase_sol_pf = solve(prob_phase_pf, BestOf(best_algs()); ranges=nodeg_ranges, kwargs...)
    # push!(phase_data_pf, phase_sol_pf)
    # repeat for level branch
    level_sol = solve(prob, BestOf(best_algs()); ranges=level_ranges, kwargs...)
    push!(level_data, level_sol)
    # repeat for level branch with penalty factor instead
    # level_sol_pf = solve(prob_level_pf, BestOf(best_algs()); ranges=nodeg_ranges, kwargs...)
    # push!(level_data_pf, level_sol_pf)
    # non-degenerate branch
    nodeg_sol = solve(prob_nodeg, BestOf(best_algs()); kwargs..., MaxTime=MaxTime / 3, initials=[2.0, ε0 + 0.1], ranges=nodeg_ranges)
    push!(nodeg_data, nodeg_sol)
end
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
    f_homogeneous = vlines!(ax, [0]; color=:red, linestyle=:dot)

    return figs, f_2s, f_homogeneous
end

branch_fig = let level_data = level_data
    with_theme(theme_latexfonts()) do
        fig = Figure(size=0.7 .* (600, 400), fontsize=15)
        g = fig[1, 1] = GridLayout()
        gl = fig[1, 2] = GridLayout()

        xticks = LinearTicks(4)
        ax1 = Axis(g[1, 1]; xlabel=L"δε_2/Δ", ylabel="LD", xticks)
        ax2 = Axis(g[2, 1]; xlabel=L"δε_2/Δ", ylabel="Excitation gap", xticks)
        linkxaxes!(ax1, ax2)
        hidexdecorations!(ax1; ticks=true, grid=false)
        CairoMakie.ylims!(ax1, (0, 1))

        fs_ld, f_ld_2s, f_ld_homogeneous = branch_plot!(ax1, δε2s, [phase_data, level_data, nodeg_data], ss2, x -> target(x.optsol))
        fs_excgap, f_excgap_2s, f_excgap_homogeneous = branch_plot!(ax2, δε2s, [phase_data, level_data, nodeg_data], ss2, x -> x.optsol.excgap)
        # fs_excgap, f_excgap_2s, f_excgap_homogeneous = branch_plot!(ax3, δε2s, [phase_data, level_data, nodeg_data], homogeneous_ss2, x -> x.optsol.gap)


        Legend(gl[1, 1],
            [fs_ld..., f_ld_2s, f_ld_homogeneous],
            [["Phase branch", "Level branch", "Topological"][1:length(fs_ld)]..., "Two-site", "Homogeneous"])
        colsize!(fig.layout, 2, 100)
        rowgap!(g, 1, 5)
        fig
    end
end
branch_fig2 = let level_data = level_data, phase_data = phase_data
    with_theme(theme_latexfonts()) do
        fig = Figure(size=0.75 .* (800, 500), fontsize=20)
        # g = fig[1:2, 1] = GridLayout()
        gl = fig[1:2, 1] = GridLayout()
        gr = fig[1:2, 2] = GridLayout()
        # grt = fig[1:2, 2] = GridLayout()
        # grb = fig[2, 2] = GridLayout()
        # gl = fig[1, 2] = GridLayout()
        xticks = LinearTicks(4)
        ax1 = Axis(gl[1, 1]; xlabel=L"δε_2/Δ", ylabel="||∇δE||")
        ax2 = Axis(gl[2, 1]; xlabel=L"δε_2/Δ", ylabel="δϕ", yticks=(pi * [0, 1 / 2, 1], ["0", L"\frac{\pi}{2}", "π"]), xticks)
        ax3 = Axis(gr[2, 1]; xlabel=L"δε_2/Δ", ylabel=L"ε_{13}", xticks)
        linkxaxes!(ax1, ax2, ax3)
        hidexdecorations!(ax1; ticks=true, grid=false)
        CairoMakie.ylims!(ax1, (0, 0.75))
        CairoMakie.ylims!(ax2, (0, pi))
        datas = [phase_data, level_data, nodeg_data]
        fs_ld, f_ld_2s, f_ld_homogeneous = branch_plot!(ax1, δε2s, datas, ss2, x -> norm(x.gradient))
        branch_plot!(ax2, δε2s, datas, ss2, x -> (x.sol[1]))
        branch_plot!(ax3, δε2s, datas, ss2, x -> (x.sol[2]))
        # branch_plot!(ax3, δε2s, [phase_data, level_data, nodeg_data], homogeneous_ss2, x -> (x.optsol.gap))

        leg = Legend(gr[1, 1],
            [fs_ld..., f_ld_2s, f_ld_homogeneous],
            ["Phase branch", "Level branch", "Topological", "Two-site sweet spot", "Homogeneous"])
        leg.tellwidth = false
        # colsize!(fig.layout, 2, 150)
        # rowsize!(gr, 2, 75)
        # colsize!(fig.layout, 2, 150)
        # rowgap!(g, 1, 5)
        fig
    end
end
branch_fig3 = let level_data = level_data, phase_data = phase_data
    with_theme(theme_latexfonts()) do
        fig = Figure(size=0.75 .* (800, 500), fontsize=20)
        gl = fig[1:2, 1] = GridLayout()
        gr = fig[1:2, 2] = GridLayout()
        xticks = LinearTicks(4)


        ax1 = Axis(gl[1, 1]; xlabel=L"δε_2/Δ", ylabel="LD", xticks)
        ax2 = Axis(gl[2, 1]; xlabel=L"δε_2/Δ", ylabel="Excitation gap", xticks)
        hidexdecorations!(ax1; ticks=true, grid=false)
        CairoMakie.ylims!(ax1, (0, 0.75))
        ax3 = Axis(gr[2, 1]; xlabel=L"δε_2/Δ", ylabel="δE/Δ", xticks)
        linkxaxes!(ax1, ax2, ax3)
        datas = [phase_data, level_data, nodeg_data]
        fs, f_2s, f_homogeneous = branch_plot!(ax1, δε2s, datas, ss2, x -> target(x.optsol))
        branch_plot!(ax2, δε2s, datas, ss2, x -> (x.optsol.excgap))
        branch_plot!(ax3, δε2s, datas, ss2, x -> (x.optsol.gap))

        leg = Legend(gr[1, 1],
            [fs..., f_2s, f_homogeneous],
            ["Phase branch", "Level branch", "Topological", "Two-site sweet spot", "Homogeneous"])
        leg.tellwidth = false
        # colsize!(fig.layout, 2, 150)
        # rowsize!(gr, 2, 75)
        # colsize!(fig.layout, 2, 150)
        # rowgap!(g, 1, 5)
        fig
    end
end
##
CairoMakie.save(plotsdir("sweet_spot_branches2.png"), branch_fig3; px_per_unit=10)
##
CairoMakie.save(plotsdir("sweet_spot_branches_appendix.png"), branch_fig2; px_per_unit=10)


##
for n in [14, 15] # Look at tuning plot for these detunings
    let k = n, hamfunc, a, data_p1, data = map(x -> x[1], nodeg_data)
        a = FermionBdGBasis(1:3)
        δε2 = δε2s[k]
        δϕ, ε1 = data[k].sol
        hamfunc = (δϕ, ε1) -> f!(cache, [δϕ, ε1, ε1 + δε2])
        εs = ε1 .+ 0.2 .* range(-1, 1, 50)
        δϕs = range(0, pi, 51)
        newdata = [(
            begin
                H = hamfunc(δϕ, ε1)
                merge(fullsolve(H, c),
                    get_gap_derivatives(x -> hamfunc(x...), c, [δϕ, ε1]))
            end
        ) for ε1 in εs, δϕ in δϕs] |> permutedims
        d = Dict("ss" => data[k], "data" => newdata, "N" => N, "labels" => ("ε1", "δϕ"), "x" => εs, "y" => δϕs)
        p1_hamfunc = (δϕ, ε1) -> LongerPoorMansMajoranas.perturbative_hamiltonian(a, 1; ε=[ε1, ε1 + δε2, ε1], δϕ=δϕ * [1, 1], t=fixedparams.t, θ=fixedparams.θ, Δ=fixedparams.Δ * [1, 1, 1], Ez=fixedparams.Ez)
        data_p1 = [(
            begin
                #H = LongerPoorMansMajoranas.perturbative_hamiltonian(a, 1; ε=[ε1, ε1 + δε2, ε1], δϕ=δϕ * [1, 1], t=fixedparams.t, θ=fixedparams.θ, Δ=fixedparams.Δ * [1, 1, 1], Ez=fixedparams.Ez)
                # H = LongerPoorMansMajoranas.perturbative_hamiltonian_homogeneous(a, 2; ε, δϕ, fixedparams...)
                merge(fullsolve(p1_hamfunc(δϕ, ε1), a),
                    get_gap_derivatives(x -> p1_hamfunc(x...), a, [δϕ, ε1]))
            end
        ) for ε1 in εs, δϕ in δϕs] |> permutedims

        data_coeffs = [(
            let
                Δ = fixedparams.Δ * [1, 1, 1]
                ε = [ε1, ε1 + δε2, ε1]
                θ = fixedparams.θ
                δϕ = δϕ * [1, 1]
                t = fixedparams.t
                Ez = fixedparams.Ez
                coeffs = stack(LongerPoorMansMajoranas.perturbative_coeffs(n; Δ, ε, θ, δϕ, t) for n in 1:2)
                Δ_aa, Δ_bb, Δ_ab, Δ_ba, t_aa, t_ab, t_ba, t_bb = collect(eachrow(coeffs))
                E1, E2, ε_ba, ε_ab, t_nn, Δ_nn = LongerPoorMansMajoranas.second_order_coeffs(1; Δ_ba, Δ_ab, t_ba, t_ab, Ez, ε, Δ)
                ε0 = [real(Ez - sqrt(Δ[n]^2 + ε[n]^2)) for n in 1:3]
                (; Δ_aa, Δ_bb, Δ_ab, Δ_ba, t_aa, t_ab, t_ba, t_bb, E1, E2, ε_ba, ε_ab, t_nn, Δ_nn, ε0)
            end
        ) for ε1 in εs, δϕ in δϕs] |> permutedims

        dp1 = Dict("ss" => data[k], "data" => data_p1, "N" => 3, "labels" => ("ε1", "δϕ"), "x" => εs, "y" => δϕs)
        p1 = plot_f(d, x -> LDbdg(x), clim_max=1, c=cgrad(:viridis, rev=true), ss_label="3-site sweet spot", legend=false)
        p2 = plot_f(dp1, x -> LDbdg(x), clim_max=1, c=cgrad(:viridis, rev=true), ss_label="3-site sweet spot", legend=:bottomright)
        # p1 = plot_f(d, x -> norm(x.gradient), clim_max=1, c=cgrad(:viridis, rev=true), ss_label="3-site sweet spot", legend=false)
        # p2 = plot_f(dp1, x -> norm(x.gradient), clim_max=1, c=cgrad(:viridis, rev=true), ss_label="3-site sweet spot", legend=:bottomright)
        # p1 = plot_f(d, x -> norm(x.gradient[2]) / LDf(x), clim_max=2, c=:redsblues, ss_label="3-site sweet spot", legend=false)
        # p2 = plot_f(dp1, x -> norm(x.gradient[2]) / LDf(x), clim_max=2, c=:redsblues, ss_label="3-site sweet spot", legend=:bottomright)
        # p1 = plot_f(d, x -> x.excgap, clim_max=.25, c=:viridis, ss_label="3-site sweet spot", legend = false)
        # p2 = plot_f(dp1, x -> x.excgap, clim_max=.25, c=:viridis, ss_label="3-site sweet spot")
        p3 = Plots.heatmap(εs, δϕs, map(x -> abs(x.t_aa[1]) - abs(x.Δ_aa[1]), data_coeffs), c=:redsblues, clims=0.1 .* (-1, 1))
        p4 = Plots.heatmap(εs, δϕs, map(x -> abs(x.Δ_nn) * angle(x.Δ_nn), data_coeffs), c=:redsblues, clims=0.001 .* (-1, 1))
        p5 = Plots.heatmap(εs, δϕs, map(x -> abs(x.t_nn) - abs(x.Δ_nn), data_coeffs), c=:redsblues, clims=0.1 .* (-1, 1))
        p6 = Plots.heatmap(εs, δϕs, map(x -> (x.ε_ba + x.ε_ab + x.ε0[2])^2 + (x.ε_ab + x.ε0[1])^2, data_coeffs), c=:redsblues, clims=0.01 .* (0, 1))


        # p = heatmap(εs, δϕs, map(x -> LDbdg(x), newdata)', c=:viridis, clims=(0, 1), xlabel="ε1", ylabel="δϕ", title="LD, 3 site, inhomogeneous")

        Plots.scatter!(p3, [ε1], [δϕ], c=:red)
        Plots.scatter!(p4, [ε1], [δϕ], c=:red)
        Plots.scatter!(p5, [ε1], [δϕ], c=:red)
        Plots.scatter!(p6, [ε1], [δϕ], c=:red)
        # p
        Plots.plot(p1, p2, p3, p4, p5, p6, plot_title="δε2 = $δε2", layout=(3, 2), size=(800, 1200)) |> display
    end
end
##

##
sweet_spots_Ez3 = []
sweet_spots_Ez2 = []
Ezs = range(1.1, 10, 10)
for Ez in Ezs
    fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0, Ez)
    f2, f2!, cache2 = hamfunc(Hδϕ_Hε(), c2, fixedparams)
    exps = range(0.1, 3, 4)
    MaxTime = 0.5
    initials = [2.2, sqrt(fixedparams.Ez^2 - fixedparams.Δ^2)]#2.95]
    homogeneous_ss2 = find_sweet_spot((f2, f2!, cache2), c2, Hδϕ_Hε(); exps, MaxTime, target, minexcgap=0, alg=BestOf(best_algs()), initials)
    f, f!, cache = hamfunc(Rδϕ_Rε(), c, fixedparams)
    fh, fh!, cacheh = hamfunc(Hδϕ_Hε(), c, fixedparams)
    homogeneous_ss = find_sweet_spot((fh, fh!, cacheh), c, Hδϕ_Hε(); exps, MaxTime, target, minexcgap=0, alg=BestOf(best_algs()), initials)
    # δϕ0 = homogeneous_ss.sol[1]
    push!(sweet_spots_Ez3, merge(homogeneous_ss, get_gap_derivatives(fh, c, homogeneous_ss.sol)))
    push!(sweet_spots_Ez2, merge(homogeneous_ss2, get_gap_derivatives(f2, c2, homogeneous_ss2.sol)))
end
##
p_Ez_excgap = plot(Ezs, [x.optsol.excgap for x in sweet_spots_Ez3], markers=true, ylabel="excgap", xlabel="Ez", ylims=(0, 0.8))
plot!(Ezs, [LDbdg(x.optsol) for x in sweet_spots_Ez2], markers=true, ls=:dash)
p_Ez_LD = plot(Ezs, [LDbdg(x.optsol) for x in sweet_spots_Ez3], markers=true, ylabel="LDbdg", xlabel="Ez", clims=(0, 1.5))
plot!(Ezs, [LDbdg(x.optsol) for x in sweet_spots_Ez2], markers=true, ls=:dash)
plot(p_Ez_excgap, p_Ez_LD)


##
inhomo_sols = []
for k in eachindex(Ezs)
    fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0, Ez=Ezs[k])
    alg = BestOf(best_algs())
    f, f!, cache = hamfunc(Rδϕ_Rε(), c, fixedparams)
    minexcgap = sweet_spots_Ez3[k].optsol.excgap * 0.1
    δϕ0, ε0 = sweet_spots_Ez3[k].sol |> collect
    prob = OptProb(; hamfunc=x -> f!(cache, [x[1], x[2], x[2] + x[3]]), basis=c, optparams=Rδϕ_Rε(), target)
    ranges = [(0.0, 1.0pi), ε0 .+ 1 .* (-1, 1), (-0.5, 0.5)]
    initials = [δϕ0, ε0, 0.0]
    sol = solve(prob, alg; minexcgap, maxiters=10000, MaxTime=0.5, exps, initials, ranges)
    push!(inhomo_sols, merge(sol, get_gap_derivatives(f, c, sol.sol)))
end
##
c_main = 1
p_Ez_excgap = plot(Ezs, [x.optsol.excgap for x in sweet_spots_Ez3], markers=true, ylabel="excgap", xlabel="Ez", ylims=(0, 0.8), c=c_main, label="3-site homogeneous")
plot!(Ezs, [LDbdg(x.optsol) for x in sweet_spots_Ez2], markers=true, label="2-site")
plot!(Ezs, [LDbdg(x.optsol) for x in inhomo_sols], markers=true, ls=:dash, c=c_main, label="3-site inhomogeneous")
p_Ez_LD = plot(Ezs, [LDbdg(x.optsol) for x in sweet_spots_Ez3], markers=true, ylabel="LDbdg", xlabel="Ez", clims=(0, 1.5), c=c_main, label="3-site homogeneous")
plot!(Ezs, [LDbdg(x.optsol) for x in sweet_spots_Ez2], markers=true, label="2-site")
plot!(Ezs, [LDbdg(x.optsol) for x in inhomo_sols], markers=true, ls=:dash, c=c_main, label="3-site inhomogeneous")
p_Ez_gap = plot(Ezs, [abs(x.optsol.gap) for x in sweet_spots_Ez3], markers=true, ylabel="gap", xlabel="Ez", clims=(1e-5, 0.1), yscale=:log, c=c_main)
plot!(Ezs, [abs(x.optsol.gap) for x in sweet_spots_Ez2], markers=true)
plot!(Ezs, [abs(x.optsol.gap) for x in inhomo_sols], markers=true, ls=:dash, c=c_main)
p_Ez_grad = plot(Ezs, [norm(x.gradient) for x in sweet_spots_Ez3], markers=true, ylabel="|gradient|", xlabel="Ez", clims=(0, 1.5), c=c_main)
plot!(Ezs, [norm(x.gradient) for x in sweet_spots_Ez2], markers=true)
plot!(Ezs, [norm(x.gradient) for x in inhomo_sols], markers=true, ls=:dash, c=c_main)
p_Ez_hess = plot(Ezs, [norm(x.hessian) for x in sweet_spots_Ez3], markers=true, ylabel="|hessian|", xlabel="Ez", clims=(0, 1.5), c=c_main)
plot!(Ezs, [norm(x.hessian) for x in sweet_spots_Ez2], markers=true)
plot!(Ezs, [norm(x.hessian) for x in inhomo_sols], markers=true, ls=:dash, c=c_main)
plot(p_Ez_excgap, p_Ez_LD, p_Ez_grad, p_Ez_hess, p_Ez_gap)
plot(p_Ez_LD, p_Ez_excgap)
##
savefig("inhomogeneous_Vz_improvement.png")

##
plot(Ezs, stack(map((h, i) -> [h.sol[1] - i.sol[1], h.sol[2] - i.sol[2], h.sol[2] - i.sol[3] - i.sol[2]], sweet_spots_Ez3, inhomo_sols))', marker=true, label=["δϕ" "ε1" "ε2"], ylabel="difference", xlabel="Ez")
##
[LDbdg(fullsolve(hamfunc(Rδϕ_Rε(), c, merge(fixedparams, (; Ez)))[1](collect(sol.sol)), c)) for (Ez, sol) in zip(Ezs, inhomo_sols)] |> plot
[LDbdg(fullsolve(hamfunc(Hδϕ_Hε(), c, merge(fixedparams, (; Ez)))[1](collect(sol.sol)[1:2]), c)) for (Ez, sol) in zip(Ezs, inhomo_sols)] |> plot!
[LDbdg(fullsolve(hamfunc(Rδϕ_Rε(), c, merge(fixedparams, (; Ez)))[1]([sol.sol[1:2]..., sol.sol[3]]), c)) for (Ez, sol) in zip(Ezs, inhomo_sols)] |> plot!


##
let sol = inhomo_sols[1]
    get_gap_derivatives(hamfunc(Rδϕ_Rε(), c, merge(fixedparams, (; Ez=3)))[1], c, [sol.sol[1:2]..., sol.sol[3]])
end