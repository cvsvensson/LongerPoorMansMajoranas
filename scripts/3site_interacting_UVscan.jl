using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra
using JLD2
using DataFrames
using LaTeXStrings
using CairoMakie
using Accessors
using ProgressMeter

##
res = (25, 25)
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
    kwargs = (; iterations=4, initials=[0.5pi, 2.9], ranges=[(0.0, 1.0pi), (2 - U - 20V, 4 + U + 20V)], MaxTime=1, local_reltol=1e-6, verbosity=0, abstol=1e-6)
    exps = range(0.0, 4.0, kwargs.iterations)
    prob_deg = ScheduledOptProb(eigfunc, target, GapPenalty(exps))
    alg = BestOf(best_algs())#BestOf(reduce(vcat, best_algs()[1:2] for k in 1:1))
    sol_nodeg = solve(prob, alg; kwargs...)
    sol_deg = solve(prob_deg, alg; kwargs...)

    return sol_deg, sol_nodeg
end
##
UViter = Iterators.product(Us, Vs)
data = @showprogress map(UV -> get_sweet_spot(UV...), UViter)
##
## Save data
wsave(datadir("final_data", "UV-tuning2.jld2"), Dict("data" => data, "Us" => Us, "Vs" => Vs, "fixedparams" => fixedparams, "N" => N, "res" => res))
## Load data
data_dict = load(datadir("final_data", "UV-tuning.jld2"));
@unpack data, Us, Vs, fixedparams, N, res = data_dict;
##
function select_degenerate_sweet_spot(data, cutoff)
    i = findfirst(x -> abs(x.optsol.gap) < cutoff, data.all_sols)
    data.all_sols[i]
end

##
cbwidth = 10
levelcolors = [:darkorange1, :crimson]
linewidth = 1.3
contourcolormap = cgrad(:managua, rev=true)
contourstyle = :dash
cbwidth = 10
linewidth = 1.3
fig_UV = let data = data
    with_theme(theme_latexfonts()) do
        fig = Figure(size=0.8 .* (600, 300), fontsize=20, figure_padding=5)

        g = fig[1, 1] = GridLayout()
        xticks = WilkinsonTicks(3)
        yticks = WilkinsonTicks(3)#(pi * [0, 1 / 2, 1], [L"0", L"\frac{\pi}{2}", L"π"])
        ax = Axis(g[1, 1]; xlabel=L"U_l/\Delta", title=L"\mathrm{mLD}_0", ylabel=L"U_\mathrm{nl}/\Delta", yticks, xticks)
        ax1 = Axis(g[1, 2]; xlabel=L"U_l/\Delta", title=L"\mathrm{mLD}", ylabel=L"U_\mathrm{nl}/\Delta", yticks, xticks)
        linkaxes!(ax, ax1)

        hideydecorations!(ax1)
        target = x -> LD_cells(x)#x -> norm(x.reduced.two_cells)

        data = map(data) do (deg, nodeg)
            (select_degenerate_sweet_spot(deg, 1e-2), nodeg)
        end

        contour_kwargs = (; color=:red, levels=[0.0], linewidth=1.3)
        LD_deg = reshape(map(ss -> LD_cells(ss[1].optsol), data), res...)
        LD_nodeg = reshape(map(ss -> LD_cells(ss[2].optsol), data), res...)
        hmap_kwargs = (; interpolate = false, colormap=Reverse(:viridis), colorscale=identity, colorrange=(0, 1.01 + 0maximum([LD_deg..., LD_nodeg...])))
        hmap = heatmap!(ax, Us, Vs, LD_deg; hmap_kwargs...)
        f_egap = heatmap!(ax1, Us, Vs, LD_nodeg; hmap_kwargs...)
        # f_Q = contour!(ax1, εs_high_res, δϕs_high_res, dataQ2; contour_kwargs...)
        # l_Q = lines!(ax1, Float64[], Float64[]; label="Q=0", contour_kwargs...)
        # axislegend(ax1; position=:lt)

        ticks = ([0, 0.25, 1 / 2, 1], ["0", ".25", "0.5", "1"])

        ticklabelsize = 16
        Colorbar(g[1, 3], hmap; width=cbwidth, ticksize=cbwidth, tickalign=true, ticks=LinearTicks(3), ticklabelsize)

        Label(g[1, 3, Bottom()], " LD", tellwidth=false, tellheight=false, fontsize=20)


        # make a), b) labels
        Label(g[1, 1, TopLeft()], "(a)", tellwidth=false, tellheight=false, fontsize=20)
        Label(g[1, 2, TopLeft()], "(b)", tellwidth=false, tellheight=false, fontsize=20)

        fig
    end
end
##
save(plotsdir("3_site_UV_sweet_spot.png"), fig_UV; px_per_unit=10)

##
heatmap(map(d -> d[1].all_sols[3].optsol.gap, data))