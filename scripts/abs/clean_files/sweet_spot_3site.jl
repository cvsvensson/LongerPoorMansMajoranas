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
# includet(scriptsdir("abs", "phase_plots.jl"))
includet(scriptsdir("abs", "phase_misc.jl"))

fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0.0, Ez=3)
## 2site sweet spot
c2 = bdg ? FermionBdGBasis(1:2, (:↑, :↓)) : FermionBasis(1:2, (:↑, :↓); qn=QuantumDots.parity)
f2, f2!, cache2 = hamfunc(Hδϕ_Hε(), c2, fixedparams)
exps = range(0.1, 4, 5)
homogeneous_ss2 = find_sweet_spot((f2, f2!, cache2), c, Hδϕ_Hε(); exps, MaxTime=1, target=LDbdg, minexcgap=0, alg=BestOf(best_algs()[1:end-1]))

##
initials = [2.2, 2.95]
bdg = true
N = 3
c = bdg ? FermionBdGBasis(1:N, (:↑, :↓)) : FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
f, f!, cache = hamfunc(Rδϕ_Rε(), c, fixedparams)
fh, fh!, cacheh = hamfunc(Hδϕ_Hε(), c, fixedparams)

exps = range(0.1, 4, 5)
homogeneous_ss = find_sweet_spot((fh, fh!, cacheh), c, Hδϕ_Hε(); exps, MaxTime=2, target=LDbdg, minexcgap=0, alg=BestOf(best_algs()[1:end-1]))
# δϕ0 = homogeneous_ss.sol[1]
ε0 = homogeneous_ss.sol[2]
initials = homogeneous_ss.sol

##
MaxTime = 0.5
minexcgap = homogeneous_ss.optsol.excgap * 0.9
data = []
ε2s = ε0 .+ 0.5 .* range(-1, 1, length=20)
for ε2 in ε2s
    alg = BestOf(best_algs()[1:end-1])
    target = bdg ? LDbdg : LD
    hamfunc = δϕε1 -> f!(cache, [δϕε1[1], δϕε1[2], ε2])
    prob = OptProb(; hamfunc, basis=c, optparams=Rδϕ_Rε(), target)
    ranges = [(0.0, 1.0pi), ε0 .+ 1 .* (-1, 1)]
    sol = solve(prob, alg; minexcgap, maxiters=100000, MaxTime, exps, initials, ranges)
    push!(data, sol)
end
##
p_LD = plot(ε2s, map(x -> LDbdg(x.optsol), data), markers=true, ylabel="LD", xlabel="ε2", label="Inhomogeneous")
scatter!([ε0], [LDbdg(homogeneous_ss.optsol)], c=:red, label="homogeneous")
hline!([LDbdg(homogeneous_ss2.optsol)], label="2 site sweet spot", lw=2, ls=:dash)
##
p_gap = plot(ε2s, map(x -> x.optsol.excgap, data), markers=true, ylabel="excgap", xlabel="ε2", label="Inhomogeneous")
scatter!([ε0], [homogeneous_ss.optsol.excgap], c=:red, label="homogeneous")
hline!([homogeneous_ss2.optsol.excgap], ls=:dash, lw=2, label="2 site sweet spot")
##
plot(p_LD, p_gap)
##
savefig("inhomogeneous_sweet_spot_3site.png")
##
fig = plot(; yscale=:log);
let
    f = d -> LDf(d.optsol)
    f = d -> LDbdgmax(d.optsol)
    ss = map(d -> d["ss"], data)
    f = d -> abs(d.optsol.gap)
    # f = d -> d.sol[1]
    ds = map(x -> f.(x["ss"].all_ss), data)
    d = map(f, ss)
    println(sum(d))
    plot!(fig, d, ls=:solid, ylims=(0.0000001, 1))
    dsp = [[(i, d) for d in ds[i]] for i in eachindex(ds)]
    println(dsp)
    scatter!(fig, vcat(dsp...))

end
fig
