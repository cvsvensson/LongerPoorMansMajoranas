using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra#, BlackBoxOptim
using Plots
using Symbolics
using Folds
using Accessors
using JLD2
using DataFrames
using LaTeXStrings
synceddir(args...) = joinpath(ENV["Dropbox"], "data", "LongerPoorMans", args...)
includet(scriptsdir("abs", "phase_plots.jl"))
includet(scriptsdir("abs", "phase_misc.jl"))
load_full_data(N, folder="high_res", dir=synceddir) = wload(dir("phase_diagram", folder, "full_N=$(N)_fixedparams=(t = 0.5, θ = QuantumDots.DiffChainParameter{Float64}(2.746801533890032), V = 0, Δ = 1, U = 0.0, Ez = 3).jld2"))
load_kitaev_data(N, folder="high_res", dir=synceddir) = wload(dir("phase_diagram", folder, "kitaev_N=$(N)_fixedparams=(t = 0.5, θ = QuantumDots.DiffChainParameter{Float64}(2.746801533890032), V = 0, Δ = 1, U = 0.0, Ez = 3).jld2"))
##
default_scale = cgrad(:viridis, rev=true, scale=x -> exp(abs(x)));
default_target = LongerPoorMansMajoranas.LDbdg
default(; c=default_scale,
    fontfamily="Computer Modern",
    colorbar_titlefontrotation=-90,
    thickness_scaling=1.5
)
##
# Ns = [2, 3, 4, 5, 20, 40]
Fdatas = Dict(zip([2:5...], load_full_data.([2:5...], "high_res", datadir)))
# Fdatas[2] = load_full_data(2, "high-res-ldbdg", datadir)
Fdatas_ss = Dict(zip(2:5, load_full_data.(2:5, "ss_bdg_noexcgap_bestof", datadir)))
foreach(N -> N in keys(Fdatas) ? (Fdatas[N]["ss"] = Fdatas_ss[N]["ss"]) : nothing, collect(keys(Fdatas_ss)))

# Fdatas_ss = Dict(zip(Ns, load_full_data.(Ns[1:4], "ss")))
# Fdatas_ss2 = Dict(zip(Ns, load_full_data.(Ns[1:4], "ss")))
# Fdatas_mpi = Dict(zip([2, 3, 4, 5, 10, 20], load_full_data.([2, 3, 4, 5, 10, 20], "mpi", datadir)))
# foreach(N -> (Fdatas[N]["ss"] = Fdatas_ss[N]["ss"]), Ns[1:4])
## Get perturbative data
a = FermionBdGBasis(1:3)
function perturbative_solutions(a, M, fixedparams, labels, x, y; transport=missing)
    xy = NamedTuple(Symbol.(labels) .=> (x, y .* ones(2)))
    fp = filter(kv -> kv[1] ∉ (:U, :V), pairs(fixedparams)) |> NamedTuple
    params = merge(fp, xy)
    H = perturbative_hamiltonian(a, M; params...)
    fullsolve(H, a; transport)
end
F3data = Fdatas[3]
εs = F3data["x"]
δϕs = F3data["y"]
labels = F3data["labels"]
fixedparams = F3data["fixedparams"]
perturbative_data = [join_data(nothing, [perturbative_solutions(a, M, fixedparams, labels, [x, x, x], y) for y in δϕs, x in εs], 3, (εs, δϕs, ("ε", "δϕ")), "perturbative", false, "") for M in 0:2];
pfdata = [perturbative_data..., F3data]

##
phase_plots = [plot_f(Fdatas[N], LDbdg; clim_max=2, leg=:bottomleft, plot_ss=!ismissing(Fdatas[N]["ss"]), title="$N sites", colorbar_title="LD", ss_label=false, c=default_scale) for N in sort(collect(keys(Fdatas)))]
phase_plots_notitle = [plot_f(Fdatas[N], LongerPoorMansMajoranas.LDbdg; leg=:bottomleft, plot_ss=!ismissing(Fdatas[N]["ss"]), title="", colorbar_title="LD", ss_label=false, c=default_scale) for N in sort(collect(keys(Fdatas)))]
# phase_plots = [plot_f(Fdatas_mpi[N], Base.Fix1(*, 1) ∘ LongerPoorMansMajoranas.LDbdgmax; c=default_scale, leg=:bottomleft, plot_ss=false, title="$N sites", colorbar_title="1-MP", ss_label=false) for N in sort(collect(keys(Fdatas_mpi)))]
##
display.(phase_plots[1:end])
##
foreach((p, N) -> savefig(p, plotsdir(string("phase_diagram_LD_$(N)", ".svg"))), phase_plots, sort(collect(keys(Fdatas))))
foreach((p, N) -> savefig(p, plotsdir(string("phase_diagram_LD_$(N)_notitle", ".svg"))), phase_plots_notitle, sort(collect(keys(Fdatas))))
##
let N = 5
    plot(plot_f(Fdatas_ss3[N], LongerPoorMansMajoranas.MPI2; leg=:bottomleft, ylabel="", colorbar=false, plot_ss=true, title="$N sites", colorbar_title="1-MP", ss_label=false),
        plot_f(Fdatas_ss3[N], MPI; leg=:bottomleft, ylabel="", plot_ss=true, colorbar=false, title="$N sites", colorbar_title="1-MP", ss_label=false))
end

##
p_full_2_40 = let
    margins = (; right_margin=0Plots.mm, left_margin=0Plots.mm)
    p2 = plot_f(Fdatas[2], default_target, leg=:bottomleft, colorbar=false, plot_ss=false, title="2 sites")
    p40 = plot_f(Fdatas[40], default_target, leg=false, ylabel="", yshowaxis=false, colorbar=false, plot_ss=false, title="40 sites")
    plot_ss!(p2, Fdatas[2]["ss"], Fdatas[2]["N"], label="2 site sweet spot")
    plot_ss!(p40, Fdatas[2]["ss"], Fdatas[2]["N"])

    h2 = scatter([0], [0]; zcolor=[0], clims=(0, 1),
        xlims=(1, 1.1), xshowaxis=false, yshowaxis=false, label="", colorbar_title="LD", grid=false)
    l = @layout [grid(1, 2) a{0.001w}]
    p_all = plot(p2, p40, h2, layout=l, size=(1000, 400))
    # plot(p2, p40, size=(1000, 400))
end
savefig(p_full_2_40, plotsdir(string("phase_diagram_LD_2_40", ".svg")))


## Perturbation theory
p_perturbation = let
    level_magnitude = 0.01
    #c = default_scale#cgrad(:viridis, rev=true, scale=x -> exp(4abs(x)))
    kwargs = (; colorbar=false, level_magnitude, fontfamily="Computer Modern",
        left_margin=-0Plots.mm, right_margin=-2Plots.mm, top_margin=2Plots.mm, bottom_margin=1Plots.mm)
    f = LongerPoorMansMajoranas.LDf
    p1 = plot_f(pfdata[2], f; legend_position=:bottomleft, plot_ss=false, title=L"H_1^\mathrm{eff}", kwargs...)
    p2 = plot_f(pfdata[3], f; legend=false, plot_ss=false, ylabel="", yshowaxis=false, title=L"H_2^\mathrm{eff}", kwargs...)
    p3 = plot_f(pfdata[4], f; legend=false, plot_ss=false, ylabel="", yshowaxis=false, title=L"H", kwargs...)
    h2 = scatter([0], [0]; zcolor=[0], clims=(0.0, 1.0),
        xlims=(1, 1.1), xshowaxis=false, yshowaxis=false, label="", colorbar_title="LD", grid=false)
    l = @layout [grid(1, 3) a{0.001w}]
    p_all = plot(p1, p2, p3, h2, layout=l, thickness_scaling=1.5, size=(900, 300))
end
savefig(p_perturbation, plotsdir(string("perturbation_theory_LDf", ".svg")))

##
## Laod Kitaev data
K2data = load_kitaev_data(2, "kitaev-high-res-ldbdg", datadir)
K40data = load_kitaev_data(40, "kitaev-high-res-ldbdg", datadir)
## Kitaev phase diagram plots
p_kitaev = let
    level_magnitude = 0.1
    c = default_scale
    kwargs = (; c, colorbar=false, level_magnitude, fontfamily="Computer Modern",
        left_margin=-0Plots.mm, right_margin=-2Plots.mm, top_margin=2Plots.mm, bottom_margin=1Plots.mm)

    f = LongerPoorMansMajoranas.LDbdg

    p1 = plot_f(K2data, f; legend_position=:bottomleft, plot_ss=false, title="2 sites", kwargs...)
    # p2 = plot_f(pfdata[3], MP; legend=false, plot_ss=false, ylabel="", yshowaxis=false, title=L"H_2^\mathrm{eff}", kwargs...)
    p3 = plot_f(K40data, f; legend=false, plot_ss=false, ylabel="", yshowaxis=false, title="40 sites", kwargs...)
    h2 = scatter([0], [0]; zcolor=[0], clims=(0, 1),
        xlims=(1, 1.1), xshowaxis=false, yshowaxis=false, label="", c, colorbar_title="LD", grid=false)
    l = @layout [grid(1, 2) a{0.001w}]
    p_all = plot(p1, p3, h2, layout=l, thickness_scaling=1.5, size=(900, 400))
end
##
savefig(p_kitaev, plotsdir(string("kitaev_phase_diagram_LD_2_40", ".svg")))

##
p_scaling = let
    Ns = sort(filter(k -> k ∉ [], collect(keys(Fdatas_ss))))
    p = plot(; xlabel="N", xticks=Ns[1]:2:Ns[end], markers=true, frame=:box, legend=:topright, size=0.9 .* (600, 400), ylims=(1e-4, 1), thickness_scaling=1.5, yscale=:log10, yticks=10.0 .^ (-4:1:1), fontfamily="Computer Modern")
    plot!(Ns, map(N -> LDbdg(Fdatas_ss[N]["ss"].optsol), Ns); c=:lighttest, legend=false, markers=true, ylabel="LD")
end
##
savefig(p_scaling, plotsdir(string("scaling_LD", ".svg")))
