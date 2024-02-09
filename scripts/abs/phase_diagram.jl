using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra#, BlackBoxOptim
using Plots
using Symbolics
using Folds
using Accessors
using JLD2
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

##
load_full_data(N) = wload(datadir("phase_diagram", "high_res", "full_N=$(N)_fixedparams=(t = 0.5, θ = QuantumDots.DiffChainParameter{Float64}(2.746801533890032), V = 0, Δ = 1, U = 0.0, Ez = 3).jld2"))
F3data = load_full_data(3)
F4data = load_full_data(4)
F5data = load_full_data(5)
F40data = load_full_data(40)
##
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
using DataFrames
synceddir(args...) = joinpath(ENV["Dropbox"], "data", "LongerPoorMans", args...)
## load all data
data = collect_results(synceddir("phase_diagram", "lengths"))
## extract sweet sweet_spots
data_gbl = groupby(data, :labels)
foreach(data -> sort!(data, :N), data_gbl)
## plot
plot(2:16, map(x -> LD(x.optsol), data_gbl[1].ss)[1:15]; xlabel="N", label="LD (stability under local perturbations)", xticks=2:2:16, markers=true, frame=:box, legend=:topright, size=0.9 .* (600, 400))
#map(x -> LD(x.optsol)^2, data_gbl[1].ss)[1:15] .|> log |> plot;
plot!(2:16, map(x -> MPU(x.optsol), data_gbl[1].ss)[1:15]; label="1 - MPU (weight of majoranas on the outermost dot)", markers=true)
plot!(2:16, map(x -> MP(x.optsol), data_gbl[1].ss)[1:15]; label="1 - MP (normalized weight of majoranas on the outermost dot)", markers=true)
##
map(x -> x.optsol.excgap, data_gbl[1].ss)[1:15] |> plot
map(x -> x.optsol.gap, data_gbl[1].ss)[1:15] |> plot


good_Ns = [n for n in data_gbl[1].N if LD(data_gbl[1][n-1, :ss].optsol) < 1 / n]
plot_MPU(data_gbl[1][9, :])
plot_MPU(data_gbl[1][10, :])
plot_MPU(data_gbl[1][18, :])

##
a = FermionBdGBasis(1:3)
function perturbative_solutions(a, M, fixedparams, labels, x, y)
    xy = NamedTuple(Symbol.(labels) .=> (x, y .* ones(2)))
    fp = filter(kv -> kv[1] ∉ (:U, :V), pairs(fixedparams)) |> NamedTuple
    params = merge(fp, xy)
    H = BdGMatrix(perturbative_hamiltonian(a, M; params...) |> Hermitian; check=false)
    fullsolve(H, a)
end
##
N = 3
fixedparams = (; t=0.001, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0.0, Ez=4)
Kdata = calculate_kitaev_phase_data(N; save=false, res=(53, 50), folder=nothing)
Fdata = calculate_full_phase_data(N; save=false, res=(53, 50), scale=25, fixedparams, optimize=false, folder=nothing)
εs = Fdata["x"]
δϕs = Fdata["y"]
##
labels = Fdata["labels"]
perturbative_data = [join_data(nothing, [perturbative_solutions(a, M, fixedparams, labels, x, y) for y in δϕs, x in εs], 3, (εs, δϕs, ("ε", "δϕ")), "perturbative", false, "") for M in 0:2];

##
plot([map(data -> plot_LD(data), perturbative_data)..., plot_LD(Fdata)]...)
##
plot([map(data -> plot_gap(data), perturbative_data)..., plot_gap(Fdata)]...)
plot([map(data -> plot_MPU(data), perturbative_data)..., plot_MPU(Fdata)]...)

