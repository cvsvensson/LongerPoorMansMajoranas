using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra
using JLD2
using DataFrames
using LaTeXStrings
using CairoMakie
using ProgressMeter
synceddir(args...) = joinpath(ENV["Dropbox"], "data", "LongerPoorMans", args...)

##
data = []
res = (500, 500)
N = 40
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0, Ez=3)
bdg = iszero(fixedparams.U) && iszero(fixedparams.V)
c = bdg ? FermionBdGBasis(1:N, (:↑, :↓)) : FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
target = LD_cells
f, f!, cache = hamfunc(Hδϕ_Hε(), c, fixedparams)
εs = sqrt(fixedparams.Ez^2 - first(fixedparams.Δ)^2) .+ range(-0.2, 0.25, length=res[1])
δϕs = range(0, pi, length=res[2])
iter = Iterators.product(δϕs, εs) |> collect
hf = δϕε -> f!(cache, δϕε)
eigfunc = δϕε -> diagonalize(f!(cache, δϕε), c)
mapfunc(δϕε) = all_info(eigfunc(δϕε))
##
data = @showprogress map(mapfunc, iter)
##
## Get sweet spots
exps = range(-1.0, 4.0, 6)
prob = ScheduledOptProb(eigfunc, target, GapPenalty(exps))
prob_level = ScheduledOptProb(eigfunc, target, GapPenalty(exps) + ScheduledPenalty((sol, x, i) -> 10^(-i)abs(get_gap(eigfunc(x .+ [0, 0.01])) - get_gap(sol)) / 0.01))
prob_phase = ScheduledOptProb(eigfunc, target, GapPenalty(exps) + ScheduledPenalty((sol, x, i) -> 10^(-i)abs(get_gap(eigfunc(x .+ [0.01, 0])) - get_gap(sol)) / 0.01))
prob_nodeg = ScheduledOptProb(eigfunc, target)
kwargs = (; iterations=length(exps), initials=[0.5pi, 2.9], ranges=[(0.0, 1.0pi), (2.5, 3.2)], MaxTime=5)
ss_deg = solve(prob, BestOf(best_algs()); kwargs...)
ss_phase = solve(prob, BestOf(best_algs()); kwargs..., ranges=[(0.0, 1.0pi), (2.7, 2.85)])
ss_level = solve(prob, BestOf(best_algs()); kwargs..., ranges=[(0.5pi, 1.0pi), (2.9, 3.1)])
ss_nodeg = solve(prob_nodeg, BestOf(best_algs()); kwargs...)
## Save data
wsave(datadir("final_data", "$N-site-tuning.jld2"), Dict("data" => data, "ss_deg" => ss_deg, "ss_nodeg" => ss_nodeg, "εs" => εs, "δϕs" => δϕs, "fixedparams" => fixedparams, "N" => N, "res" => res, "target" => target, "bdg" => bdg))
