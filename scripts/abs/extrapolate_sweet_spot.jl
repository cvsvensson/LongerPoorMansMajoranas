# In this file, we extrapolate the sweet spot for the N=2 case to larger chains. If Δ is large, this should still be the sweet spot
using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra
using Plots
using Symbolics
using Folds
using Accessors

function extrapolate_ss_initials(sol, N)
    δϕ = sol[1]
    vcat([δϕ for k in 1:div(N, 2)], [sol[2] for n in 1:div(N + 1, 2)])
end
function extrapolate_ss(sol, fixedparams, N)
    δϕ = sol[1]
    ϕ = δϕ .* (0:N-1)
    ε = parameter(sol[2], :homogeneous)
    Δ = @. exp(1im * ϕ) * fixedparams.Δ
    newparams = merge(fixedparams, (; Δ, ε))
    basis = FermionBdGBasis(1:N, (:↑, :↓))
    h = LongerPoorMansMajoranas.whamiltonian(basis; newparams...) |> hermitianpart! |> BdGMatrix
    fullsolve(h, basis)
end

## Find sweet spot in the N=2 case
Δ0 = 1
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=Δ0, U=0.0)
N = 2
basis = FermionBdGBasis(1:N, (:↑, :↓))
optim = Rδϕ_Rε()
alginds = [1, 1, 2, 3, 4]
Ezs = collect(range(0.1, 8, length=50)) .* Δ0
Ezsols = let
    Ezsols = Vector{Any}(undef, length(Ezs))
    for (k, Ez) in collect(enumerate(Ezs))
        fp = merge(fixedparams, (; Ez))
        f, f!, cache = hamfunc(optim, basis, fp)
        prob = OptProb(; hamfunc=x -> f!(cache, x), basis, optparams=optim)
        algsols = [solve(prob, alg; minexcgap=0.1, maxiters=1000000, MaxTime=2, exps=range(0.1, 3, 5)) for alg in best_algs()[alginds]]
        Ezsols[k] = algsols
    end
    Ezsols
end;
bestsols = map(x -> findmin(y -> MPU((y.optsol)), x)[2], Ezsols)
result = Dict("optparams" => map((s, n) -> collect(s[n].sol), Ezsols, bestsols), "Ez" => Ezs, "fixedparams" => fixedparams)
##
# result2 = result
N = 3
basis = FermionBdGBasis(1:N, (:↑, :↓))
Ezs = result2["Ez"]
Ezsols3 = let
    Ezsols = Vector{Any}(undef, length(Ezs))
    for (k, Ez) in collect(enumerate(Ezs))
        fp = merge(result2["fixedparams"], (; Ez))
        f, f!, cache = hamfunc(optim, basis, fp)
        prob = OptProb(; hamfunc=x -> f!(cache, x), basis, optparams=optim)
        initials = extrapolate_ss_initials(result2["optparams"][k], N)
        println(initials)
        algsols = [solve(prob, alg; initials, minexcgap=0.1, maxiters=1000000, MaxTime=2, exps=range(0.1, 3, 5)) for alg in best_algs()[alginds]]
        Ezsols[k] = algsols
    end
    Ezsols
end;
bestsols3 = map(x -> findmin(y -> MPU((y.optsol)), x)[2], Ezsols3)
result3 = Dict("optparams" => map((s, n) -> collect(s[n].sol), Ezsols3, bestsols3), "Ez" => Ezs, "fixedparams" => fixedparams)
##
map(MPU, result3)

##
wsave(sdatadir("Ezsols_small_N3.jld2"), result3)
##
data_small2 = wload(datadir("Ezsols_small.jld2"))
data_small3 = wload(datadir("Ezsols_small_N3.jld2"))
data_big = wload(datadir("Ezsols_big3.jld2"))
##
optim = Rδϕ_Rε()
sols3 = let data = data_small3, basis = FermionBdGBasis(1:3, (:↑, :↓))
    Ezs = data["Ez"]
    fixedparams = data["fixedparams"]
    optparams = data["optparams"]
    sols = Vector{Any}(undef, length(Ezs))
    for (k, (Ez, ps)) in collect(enumerate(zip(Ezs, optparams)))
        fp = merge(fixedparams, (; Ez))
        f, f!, cache = hamfunc(optim, basis, fp)
        H = f!(cache, ps)
        sols[k] = fullsolve(H, basis)
    end
    sols
end;
##
plot(map(MPU, sols2))
plot!(map(MPU, sols3))
plot(map(LD, sols2))
plot!(map(LD, sols3))
plot(map(x->x.excgap, sols2))
plot!(map(x->x.excgap, sols3))
plot(map(x->x.gap, sols2))
plot!(map(x->x.gap, sols3))

##

## Study the extrapolated sweet spot
function get_extrapolated_sols(optdata, Ns)
    [extrapolate_ss(optdata["optparams"][k], merge(optdata["fixedparams"], (; Ez)), N) for N in Ns, (k, Ez) in enumerate(optdata["Ez"])]
end
##
Ns = 2:40
sols_small = get_extrapolated_sols(data_small, Ns)
sols_big = get_extrapolated_sols(data_big, Ns)

## Plotting
function plot_extrapolated_ss(optdata, extrapolated_ss; clims=(-3, 0))
    Ezs = optdata["Ez"]
    sols = extrapolated_ss
    c = cgrad(:viridis, rev=true)
    xlabel = "Ez/Δ"
    ylabel = "N"
    pmpu = heatmap(Ezs, Ns, map(log10 ∘ MPU, sols); c, clims, xlabel, ylabel, title="log10(1-MPU)")
    pld = heatmap(Ezs, Ns, map(log10 ∘ LD, sols); c, clims, xlabel, ylabel, title="log10(LD)")
    pgap = heatmap(Ezs, Ns, map(x -> log10(abs(x.gap)), sols); c, clims, xlabel, ylabel, title="log10(Egap)")
    pexcgap = heatmap(Ezs, Ns, map(x -> (abs(x.excgap)), sols); c, clims=(0, 0.2), xlabel, ylabel, title="Excgap")
    plot(pmpu, pld, pgap, pexcgap; plot_title="Δ = $(round(optdata["fixedparams"].Δ/optdata["fixedparams"].t, digits = 3))t", size=(800, 600))
end
##
plot_extrapolated_ss(data_small, sols_small)
plot_extrapolated_ss(data_big, sols_big)
##
clims = (-3, 0)
ps = let data = data_small, sols = sols_small
    Ezs = data["Ez"]
    pmpu = heatmap(Ezs, Ns, map(log10 ∘ MPU, sols), c=cgrad(:viridis, rev=true), clims=clims, xlabel="Ez/Δ", ylabel="N", title="log10(1-MPU)")
    pld = heatmap(Ezs, Ns, map(log10 ∘ LD, sols), c=cgrad(:viridis, rev=true), clims=clims, xlabel="Ez/Δ", ylabel="N", title="log10(LD)")
end
pmpu = heatmap(Ezs_big ./ Δ0, Ns, map(log10 ∘ abs ∘ MPU, sols_big), c=cgrad(:viridis, rev=true), clims=clims, xlabel="Ez/Δ", ylabel="N", title="log10(1-MPU)")
##
pld = heatmap(Ezs, Ns, map(log10 ∘ LD, sols), c=cgrad(:viridis, rev=true), clims=clims, xlabel="Ez/Δ", ylabel="N", title="log10(LD)")
pld = heatmap(Ezs_big ./ Δ0, Ns, map(log10 ∘ LDf, sols_big), c=cgrad(:viridis, rev=true), clims=clims, xlabel="Ez/Δ", ylabel="N", title="log10(LD)")
##
pgap = heatmap(Ezs, Ns, map(x -> log10(abs(x.gap)), sols), c=cgrad(:viridis, rev=true), clims=(-5, 0), xlabel="Ez/Δ", ylabel="N", title="log10(Egap)")
pgap = heatmap(Ezs_big, Ns, map(x -> log10(abs(x.gap)), sols_big), c=cgrad(:viridis, rev=true), clims=clims, xlabel="Ez/Δ", ylabel="N", title="log10(Egap)")
##
pexcgap = heatmap(Ezs, Ns, map(x -> (abs(x.excgap)), sols), c=cgrad(:viridis, rev=true), clims=(0, 0.2), xlabel="Ez/Δ", ylabel="N", title="Excgap")
pexcgap = heatmap(Ezs_big, Ns, map(x -> (abs(x.excgap)), sols_big), c=cgrad(:viridis, rev=true), clims=(0, 0.2), xlabel="Ez/Δ", ylabel="N", title="Excgap")
##
heatmap(Ezs, Ns, map(x -> log10(abs(x.gapratio)), sols_big), c=cgrad(:viridis, rev=true), clims=(-20, 0), xlabel="Ez/Δ", ylabel="N", title="gapratio")
##
plot(pmpu, pld, pgap, pexcgap; plot_title="Δ = 20t", size=(800, 600))