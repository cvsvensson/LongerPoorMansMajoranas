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

fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0.0, Ez=1.5)
## 2site sweet spot
bdg = true
c2 = bdg ? FermionBdGBasis(1:2, (:↑, :↓)) : FermionBasis(1:2, (:↑, :↓); qn=QuantumDots.parity)
f2, f2!, cache2 = hamfunc(Hδϕ_Hε(), c2, fixedparams)
exps = range(0.1, 4, 5)
homogeneous_ss2 = find_sweet_spot((f2, f2!, cache2), c2, Hδϕ_Hε(); exps, MaxTime=2, target=LDbdg, minexcgap=0, alg=BestOf(best_algs()[1:end-1]))

##
initials = [2.2, 2.95]
N = 3
c = bdg ? FermionBdGBasis(1:N, (:↑, :↓)) : FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
f, f!, cache = hamfunc(Rδϕ_Rε(), c, fixedparams)
fh, fh!, cacheh = hamfunc(Hδϕ_Hε(), c, fixedparams)

f([1,1,1]) - fh([1,1]) |> norm

exps = range(0.1, 4, 5)
homogeneous_ss = find_sweet_spot((fh, fh!, cacheh), c, Hδϕ_Hε(); exps, MaxTime=2, target=LDbdg, minexcgap=0, alg=BestOf(best_algs()[1:end-1]))
# δϕ0 = homogeneous_ss.sol[1]
ε0 = homogeneous_ss.sol[2]
initials = homogeneous_ss.sol

##
MaxTime = 2
minexcgap = homogeneous_ss.optsol.excgap * 0.5
data = []
ε2s = ε0 .+ 0.3 .* range(-.5, 1, length=10)
for ε2 in ε2s
    alg = BestOf(best_algs()[1:end-1])
    target = bdg ? LDbdg : LD
    hamfunc = δϕε1 -> f!(cache, [δϕε1[1], δϕε1[2], ε2])
    prob = OptProb(; hamfunc, basis=c, optparams=Rδϕ_Rε(), target)
    ranges = [(0.0, 1.0pi), ε0 .+ 1 .* (-1, 1)]
    sol = solve(prob, alg; minexcgap, maxiters=10000, MaxTime, exps, initials, ranges)
    push!(data, sol)
end
##
p_LD = plot(ε2s, map(x -> LDbdg(x.optsol), data), markers=true, ylabel="LD", xlabel="ε2", label="Inhomogeneous")
scatter!([ε0], [LDbdg(homogeneous_ss.optsol)], c=:red, label="homogeneous")
hline!([LDbdg(homogeneous_ss2.optsol)], label="2 site sweet spot", lw=2, ls=:dash)
##
p_excgap = plot(ε2s, map(x -> x.optsol.excgap, data), markers=true, ylabel="excgap", xlabel="ε2", label="Inhomogeneous", clims=(0, 0.3))
scatter!([ε0], [homogeneous_ss.optsol.excgap], c=:red, label="homogeneous")
hline!([homogeneous_ss2.optsol.excgap], ls=:dash, lw=2, label="2 site sweet spot")
##
p_gap = plot(ε2s, map(x -> x.optsol.gap, data), markers=true, ylabel="gap", xlabel="ε2", label="Inhomogeneous")
scatter!([ε0], [homogeneous_ss.optsol.gap], c=:red, label="homogeneous")
hline!([homogeneous_ss2.optsol.gap], ls=:dash, lw=2, label="2 site sweet spot")
##
plot(p_LD, p_excgap, p_gap)
##
savefig("inhomogeneous_sweet_spot_3site.png")
##
include(scriptsdir("abs", "phase_plots.jl"))
let k = 7, hamfunc, a, data_p1
    a = FermionBdGBasis(1:3)
    ε2 = ε2s[k]
    δϕ, ε1 = data[k].sol
    hamfunc = (δϕ, ε1) -> f!(cache, [δϕ, ε1, ε2])
    εs = ε1 .+ 0.2 .* range(-1, 1, 100)
    δϕs = range(0, pi, 101)
    newdata = [(
        begin
            H = hamfunc(δϕ, ε1)
            fullsolve(H, c)
        end
    ) for ε1 in εs, δϕ in δϕs] |> permutedims
    d = Dict("ss" => data[k], "data" => newdata, "N" => N, "labels" => ("ε1", "δϕ"), "x" => εs, "y" => δϕs)
    data_p1 = [(
        begin
            H = LongerPoorMansMajoranas.perturbative_hamiltonian(a, 1; ε=[ε1, ε2, ε1], δϕ=δϕ * [1, 1], t=fixedparams.t, θ=fixedparams.θ, Δ=fixedparams.Δ * [1, 1, 1], Ez=fixedparams.Ez)
            # H = LongerPoorMansMajoranas.perturbative_hamiltonian_homogeneous(a, 2; ε, δϕ, fixedparams...)
            fullsolve(H, a)
        end
    ) for ε1 in εs, δϕ in δϕs] |> permutedims
    dp1 = Dict("ss" => data[k], "data" => data_p1, "N" => 3, "labels" => ("ε1", "δϕ"), "x" => εs, "y" => δϕs)
    p1 = plot_f(d, x -> LDbdg(x), clim_max=1, c=:viridis, ss_label="3-site sweet spot", legend = false)
    p2 = plot_f(dp1, x -> LDbdg(x), clim_max=1, c=:viridis, ss_label="3-site sweet spot")
    # p1 = plot_f(d, x -> x.excgap, clim_max=.25, c=:viridis, ss_label="3-site sweet spot", legend = false)
    # p2 = plot_f(dp1, x -> x.excgap, clim_max=.25, c=:viridis, ss_label="3-site sweet spot")
    # p = heatmap(εs, δϕs, map(x -> LDbdg(x), newdata)', c=:viridis, clims=(0, 1), xlabel="ε1", ylabel="δϕ", title="LD, 3 site, inhomogeneous")

    # scatter!(p, [ε1], [δϕ], c=:red)
    # p
    plot(p1, p2)
end

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
