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
##
K2Data = wload(datadir("phase_diagram", "kitaev_N=2_fixedparams=(Δ = 1,).jld2"))
K40Data = wload(datadir("phase_diagram", "kitaev_N=40_fixedparams=(Δ = 1,).jld2"))
plot_data(K2Data)
plot_data(K40Data)

##
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0.0, Ez=3)
@time data = calculate_full_phase_data(40; save=true, res=(100, 100), fixedparams, MaxTime=20, optimize=true, exps=range(0.1, 4, 5), folder="high_res")
@time Kdata = calculate_kitaev_phase_data(2; save=false, res=(50, 50))
##
for N in 2:60
    calculate_kitaev_phase_data(N; save=true, res=(50, 50))
    calculate_full_phase_data(N; save=true, res=(50, 50), fixedparams, MaxTime=10, optimize=true, exps=range(0.1, 3, 5))
end
##
plot_LD(data)
plot_LD(Kdata)
## Full model
calculate_full_phase_data(40; save=true, res=(250, 250), fixedparams, MaxTime=2, optimize=true, exps=range(0.1, 3, 5), scale=0.75, folder="high_res")
##
load_full_data(N) = wload(datadir("phase_diagram", "high_res", "full_N=$(N)_fixedparams=(t = 0.5, θ = QuantumDots.DiffChainParameter{Float64}(2.746801533890032), V = 0, Δ = 1, U = 0.0, Ez = 3).jld2"))
F2data = load_full_data(2)
F3data = load_full_data(3)
F4data = load_full_data(4)
F5data = load_full_data(5)
F40data = load_full_data(40)
F20data = load_full_data(20)
##
a = FermionBasis(1:3, qn=QuantumDots.parity)
function perturbative_solutions(a, M, fixedparams, labels, x, y; transport=missing)
    xy = NamedTuple(Symbol.(labels) .=> (x, y .* ones(2)))
    fp = filter(kv -> kv[1] ∉ (:U, :V), pairs(fixedparams)) |> NamedTuple
    params = merge(fp, xy)
    H = perturbative_hamiltonian(a, M; params...)
    fullsolve(H, a; transport)
end
εs = F3data["x"]
δϕs = F3data["y"]
labels = F3data["labels"]
perturbative_data = [join_data(nothing, [perturbative_solutions(a, M, fixedparams, labels, [x, x, x], y) for y in δϕs, x in εs], 3, (εs, δϕs, ("ε", "δϕ")), "perturbative", false, "") for M in 0:2];
pfdata = [perturbative_data..., F3data]
##
plot_f(F2data, MP; c=cgrad(:viridis, rev=true, scale=x -> exp(5x)), leg=:bottomleft, plot_ss=false, fontfamily="Computer Modern", title="2", colorbar_titlefontrotation=-90, colorbar_title="1-MP")
plot_f(F4data, MP; c=cgrad(:viridis, rev=true, scale=x -> exp(4x)), plot_ss=false, fontfamily="Computer Modern", title="4", colorbar_titlefontrotation=-90, leg=:bottomleft, colorbar_title="1-MP")
plot_f(F5data, MP; c=cgrad(:viridis, rev=true, scale=x -> exp(4x)), plot_ss=false, fontfamily="Computer Modern", title="5", colorbar_titlefontrotation=-90, leg=:bottomleft, colorbar_title="1-MP")
plot_f(F20data, MP; c=cgrad(:viridis, rev=true, scale=x -> exp(4x)), plot_ss=false, fontfamily="Computer Modern", title="20", colorbar_titlefontrotation=-90, leg=:bottomleft, colorbar_title="1-MP")
plot_f(F40data, MP; c=cgrad(:viridis, rev=true, scale=x -> exp(6x)), plot_ss=false, fontfamily="Computer Modern", title="40", colorbar_titlefontrotation=-90, leg=:bottomleft, colorbar_title="1-MP")

plot_f(F3data, MP; c=cgrad(:viridis, rev=true, scale=x -> exp(4x)), plot_ss=false, fontfamily="Computer Modern", title="Three sites", colorbar_titlefontrotation=-90, colorbar_title="1-MP")
##
let
    level_magnitude = 0.01
    c = cgrad(:viridis, rev=true, scale=x -> exp(4abs(x)))
    kwargs = (; c, colorbar=false, level_magnitude, fontfamily="Computer Modern",
        left_margin=-0Plots.mm, right_margin=-2Plots.mm, top_margin=2Plots.mm, bottom_margin=1Plots.mm)
    p1 = plot_f(pfdata[2], MP; legend_position=:bottomleft, plot_ss=false, title=L"H_1^\mathrm{eff}", kwargs...)
    p2 = plot_f(pfdata[3], MP; legend=false, plot_ss=false, ylabel="", yshowaxis=false, title=L"H_2^\mathrm{eff}", kwargs...)
    p3 = plot_f(pfdata[4], MP; legend=false, plot_ss=false, ylabel="", yshowaxis=false, title=L"H", kwargs...)
    h2 = scatter([0], [0]; zcolor=[0], clims=(0.0, 1.0),
        xlims=(1, 1.1), xshowaxis=false, yshowaxis=false, label="", c, colorbar_title="1-MP", grid=false)
    l = @layout [grid(1, 3) a{0.001w}]
    p_all = plot(p1, p2, p3, h2, layout=l, thickness_scaling=1.5, size=(900, 300))
end

##
plot_f(pfdata[4], MP; fontfamily="Computer Modern", leg=:bottomleft, plot_ss=false, colorbar_titlefontrotation=-90, size=0.6 .* (800, 600), colorbar_title="1-MP", c=cgrad(:viridis, rev=true, scale=x -> exp(4abs(x))))

##
plot_f(F3data, LD)
plot_f(F3data, MPU)
plot_LD(F3data)
plot_MPU(F4data)
plot_MPU(F5data)
##
plot_LD(F40data)
plot_MPU(F40data)
plot_MP(F40data)
plot_gap(F40data)
##
F40data_new_ss = deepcopy(F40data);
##
ss = find_sweet_spot(40; MaxTime=30, exps=range(0.1, 4, 5))
##
F40data_new_ss["ss"] = ss
plot_LD(F40data_new_ss)
##

## load all data
data = collect_results(synceddir("phase_diagram", "lengths"))
## extract sweet sweet_spots
data_gbl = groupby(data, :labels)
foreach(data -> sort!(data, :N), data_gbl)
## plot
p = plot(; xlabel="N", xticks=2:2:16, markers=true, frame=:box, legend=:topright, size=0.9 .* (600, 400), thickness_scaling=1.3, yscale=:log10, yticks=10.0 .^ (-6:2:0), fontfamily="Computer Modern");
plot!(2:16, map(x -> MP(x.optsol), data_gbl[1].ss)[1:15]; legend=false, markers=true, ylabel="1-MP")
plot(2:16, map(x -> LD(x.optsol), data_gbl[1].ss)[1:15]; label="LD (stability under local perturbations)", markers=true)
#map(x -> LD(x.optsol)^2, data_gbl[1].ss)[1:15] .|> log |> plot;
plot!(2:16, map(x -> MPU(x.optsol), data_gbl[1].ss)[1:15]; label="1 - MPU (weight of majoranas on the outermost dot)", markers=true)
##
map(x -> x.optsol.excgap, data_gbl[1].ss)[1:15] |> plot
map(x -> x.optsol.gap, data_gbl[1].ss)[1:15] |> plot


good_Ns = [n for n in data_gbl[1].N if LD(data_gbl[1][n-1, :ss].optsol) < 1 / n]
plot_MPU(data_gbl[1][9, :])
plot_MPU(data_gbl[1][10, :])
plot_MPU(data_gbl[1][18, :])

##
# a = FermionBdGBasis(1:3)
a = FermionBasis(1:3, qn=QuantumDots.parity)
function perturbative_solutions(a, M, fixedparams, labels, x, y; transport=missing)
    xy = NamedTuple(Symbol.(labels) .=> (x, y .* ones(2)))
    fp = filter(kv -> kv[1] ∉ (:U, :V), pairs(fixedparams)) |> NamedTuple
    params = merge(fp, xy)
    # H = BdGMatrix(perturbative_hamiltonian(a, M; params...) |> Hermitian; check=false)
    H = perturbative_hamiltonian(a, M; params...)
    fullsolve(H, a; transport)
end
##
N = 3
# fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0.0, Ez=3)
# fixedparams = (; t=0.879, θ=parameter(2atan(0.995962), :diff), V=0, Δ=1, U=0.0, Ez=1.3143)
transport = Transport(QuantumDots.PauliSystem, (; T=1 / 20, μ=(0.0, 0.0)))
Kdata = calculate_kitaev_phase_data(N; save=false, res=(53, 50), folder=nothing)
Fdata = calculate_full_phase_data(N; save=false, res=(53, 50), scale=1, fixedparams, optimize=true, folder=nothing, transport)
εs = Fdata["x"]
δϕs = Fdata["y"]
##
labels = Fdata["labels"]
perturbative_data = [join_data(nothing, [perturbative_solutions(a, M, fixedparams, labels, [x, x, x], y; transport) for y in δϕs, x in εs], 3, (εs, δϕs, ("ε", "δϕ")), "perturbative", false, "") for M in 0:2];

pfdata = [perturbative_data..., Fdata]
##
# plot(map(data -> plot_f(data, MPU), pfdata)[2:4]..., layout=(1, 3))
using LaTeXStrings
p1 = plot(plot_f(pfdata[2], MPU), colorbar=false, yshowaxis=true, title=L"H_1")
p2 = plot(plot_f(pfdata[3], MPU), colorbar=false, yshowaxis=false, title=L"H_2", ylabel="")
p3 = plot(plot_f(pfdata[4], MPU), colorbar=false, yshowaxis=false, title=L"H_{full}", ylabel="")

h2 = scatter([0, 0], [0, 1], zcolor=[0, 3], clims=(0, 1),
    xlims=(1, 1.1), xshowaxis=false, yshowaxis=false, label="", c=:viridis, colorbar_title="MPU", grid=false)
l = @layout [grid(1, 3) a{0.01w}]
p_all = plot(p1, p2, p3, h2, layout=l, size=(800, 300))
##
plot(map(data -> plot_f(data, sqrt ∘ MP), pfdata)...)
plot(map(data -> plot_f(data, MPU), pfdata)...)
plot(map(data -> plot_f(data, LDf), pfdata)...)
plot(map(data -> plot_f(data, x -> norm(x.reduced.cells2)), pfdata)...)
plot(map(data -> plot_f(data, x -> x.conductance[1, 1]; clims=(-0.1, 30)), pfdata)...)
plot(map(data -> plot_f(data, x -> x.conductance[1, 2]; clims=(-5, 5)), pfdata)...)
##
c = FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
T = 1 / 100
Vs = range(-0.3, 0.3, 42)
cs = conductance_sweep(c, fixedparams, F3data["ss"].sol, εs[1:6:end], Vs, T)
csp = pert_conductance_sweep(fixedparams, F3data["ss"].sol, εs[1:6:end], -Vs, T)
cs2 = conductance_sweep(c, fixedparams, F3data["ss"].sol .+ [-0.5, -0.15], εs[1:6:end], Vs, T)
csp2 = pert_conductance_sweep(fixedparams, F3data["ss"].sol .+ [-0.5, -0.15], εs[1:6:end], -Vs, T)
##
get_cs_heatmaps(data; kwargs...) = [heatmap(map(x -> x.conductance[1, 1], data.twol |> permutedims); xlabel="ε12", ylabel="Vl", kwargs...),
    heatmap(map(x -> x.conductance[2, 2], data.twor |> permutedims); xlabel="ε12", ylabel="Vr", kwargs...),
    heatmap(map(x -> x.conductance[1, 1], data.threel |> permutedims); xlabel="ε123", ylabel="Vl", kwargs...)]
hs = stack(get_cs_heatmaps.([csp..., cs], c=cgrad(:amp, scale=:exp), colorbar=false, ticks=false, clims=(0, 40)))
plot(hs..., size=(800, 800))
##
gaps = [map(pdata -> map(x -> log(abs(x.gap)), pdata["data"][1200:1400]), pfdata)...]

plot(stack((vec.(gaps))), ylims=(-18, -7))

##
N = 3
ts = exp.(range(-10, 1, 10))
fixedparams_t = [(; t, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0.0, Ez=4) for t in ts]
Fdata_t = [calculate_full_phase_data(N; save=false, res=(5, 4), scale=1 / fixedparams.t, fixedparams, optimize=false, folder=nothing) for fixedparams in fixedparams_t]
εs = Fdata_t[1]["x"]
δϕs = Fdata_t[1]["y"]

##
labels = Fdata_t[1]["labels"]
perturbative_data_t = [[join_data(nothing, [perturbative_solutions(a, M, fixedparams, labels, x, y) for y in δϕs, x in εs], 3, (εs, δϕs, ("ε", "δϕ")), "perturbative", false, "") for M in 0:2] for fixedparams in fixedparams_t];
##
data_t = [[perturbative_data..., Fdata] for (perturbative_data, Fdata) in zip(perturbative_data_t, Fdata_t)]
##
gaps = data_t .|> x -> x .|> x -> x["data"][:, 1] .|> x -> x.gap
# display.(plot.(stack.(gaps)))
##
scalings = map(g -> map(gp -> log10(norm(gp .- g[4])), g[1:3]), (gaps)) |> stack |> permutedims
plot(log10.(ts), scalings) |> display
diff(scalings; dims=1) / diff(log10.(ts))[1] |> plot