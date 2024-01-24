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
@time data = calculate_full_phase_data(4; save=false, res=(50, 50), fixedparams, MaxTime=1, optimize=true, exps=range(0.1, 3, 5))
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
F2Data = wload(datadir("phase_diagram", "full_N=2_fixedparams=(t = 0.5, θ = QuantumDots.DiffChainParameter{Float64}(2.746801533890032), V = 0, Δ = 1, U = 0.0, Ez = 3).jld2"))
F40Data = wload(datadir("phase_diagram", "full_N=40_fixedparams=(t = 0.5, θ = QuantumDots.DiffChainParameter{Float64}(2.746801533890032), V = 0, Δ = 1, U = 0.0, Ez = 3).jld2"))
F3Data = wload(datadir("phase_diagram", "full_N=3_fixedparams=(t = 0.5, θ = QuantumDots.DiffChainParameter{Float64}(2.746801533890032), V = 0, Δ = 1, U = 0.0, Ez = 3).jld2"))
##
plot_LD(F2Data)
plot_MPU(F2Data)
plot_MP(F2Data)
plot_gap(F2Data)
##
plot_LD(F40Data)
plot_MPU(F40Data)
plot_MP(F40Data)
plot_gap(F40Data)
##

using DataFrames
synceddir(args...) = joinpath(ENV["Dropbox"], "data", "LongerPoorMans", args...)
## load all data
data = collect_results(synceddir("phase_diagam", "length"))
## extract sweet sweet_spots


## plot LD and MP for sweet spots
