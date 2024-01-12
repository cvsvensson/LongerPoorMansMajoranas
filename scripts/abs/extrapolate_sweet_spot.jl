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
##
sols = Dict()
bestsols = Dict()
## Find sweet spot in the N=2 case
Δ0 = 1
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=Δ0, U=0.0)
N = 2
basis = FermionBdGBasis(1:N, (:↑, :↓))
optim = Rδϕ_Rε()
alginds = [1, 1, 2, 3, 4]
Ezs = collect(range(0.1, 7, length=30)) .* Δ0
let
    Ezsols = Vector{Any}(undef, length(Ezs))
    for (k, Ez) in collect(enumerate(Ezs))
        fp = merge(fixedparams, (; Ez))
        f, f!, cache = hamfunc(optim, basis, fp)
        prob = OptProb(; hamfunc=x -> f!(cache, x), basis, optparams=optim)
        algsols = [solve(prob, alg; minexcgap=0.1, maxiters=1000000, MaxTime=2, exps=range(0.1, 5, 5)) for alg in best_algs()[alginds]]
        Ezsols[k] = algsols
    end
    sols[2] = Ezsols
end;
##
bestsols[2] = [sol[findmin(y -> MPU(y.optsol), sol)[2]] for sol in sols[2]]
result = Dict("optparams" => map((s, n) -> collect(s[n].sol), Ezsols, bestsols), "Ez" => Ezs, "fixedparams" => fixedparams)
##
p = plot(map(x -> x[1].optsol.gap, sols[2]))
for k in Iterators.drop(eachindex(alginds), 1)
    plot!(p, map(x -> x[k].optsol.gap, sols[2]))
end
p
##
p = plot(map(x -> MPU(x[1].optsol), sols[2]))
for k in Iterators.drop(eachindex(alginds), 1)
    plot!(p, map(x -> MPU(x[k].optsol), sols[2]))
end
p
##
extrapolatedsols = [extrapolate_ss(bestsols[2][k].sol, merge(fixedparams, (; Ez)), 3) for (k, Ez) in enumerate(Ezs)]
##
N = 3
basis = FermionBdGBasis(1:N, (:↑, :↓))
let
    Ezsols = Vector{Any}(undef, length(Ezs))
    for (k, Ez) in collect(enumerate(Ezs))
        sol = bestsols[2][k]
        exsol = extrapolatedsols[k]
        fp = merge(fixedparams, (; Ez))
        f, f!, cache = hamfunc(optim, basis, fp)
        prob = OptProb(; hamfunc=x -> f!(cache, x), basis, optparams=optim)
        initials = extrapolate_ss_initials(sol.sol, N)
        println(initials)
        algsols = [solve(prob, alg; initials, minexcgap=0.05, maxiters=1000000, MaxTime=10, exps=range(0.1, 5, 5)) for alg in best_algs()[alginds]]
        Ezsols[k] = algsols
    end
    sols[3] = Ezsols
end;
bestsols[3] = [sol[findmin(y -> MPU(y.optsol), sol)[2]] for sol in sols[3]]
# result3 = Dict("optparams" => map((s, n) -> collect(s[n].sol), Ezsols3, bestsols3), "Ez" => Ezs, "fixedparams" => fixedparams)
##

##
plot(map(x -> MPU(x.optsol), bestsols[2]))
plot!(map(x -> MPU(x.optsol), bestsols[3]))
plot!(map(MPU, extrapolatedsols))
##
plot(map(x -> LD(x.optsol), bestsols[2]))
plot!(map(x -> LD(x.optsol), bestsols[3]))
plot!(map(LD, extrapolatedsols))
##
plot(map(x -> x.optsol.excgap, bestsols[2]))
plot!(map(x -> x.optsol.excgap, bestsols[3]))
plot!(map(x -> x.excgap, extrapolatedsols))
##
plot(map(x -> x.optsol.gap, bestsols[2]))
plot!(map(x -> x.optsol.gap, bestsols[3]))
plot!(map(x -> x.gap, extrapolatedsols))
##
map(MPU, result3)
##
Δ0 = 10
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=Δ0, U=0.0)
f, f!, cache = hamfunc(optim, basis, merge(fixedparams, (; Ez=5 * Δ0)))
##
ε1 = range(4.88, 4.92, length=100) .* Δ0
ε2 = range(4.85, 5, length=200) .* Δ0
δϕ = range(pi/4, π/2, length=5)
sols3d = [fullsolve(f!(cache, [δϕ, ε1, ε2]), basis) for (δϕ, ε1, ε2) in Iterators.product(δϕ, ε1, ε2)]
##
function slice_plot(sols, ε1=ε1, ε2=ε2)
    cs = [:black, :white, :black]
    ls = ([[-1e-2], [0.0], [1e-2]])
    contours = map(x -> x.gap, sols |> permutedims)

    z = map(log10 ∘ MPU, sols |> permutedims)
    clims = (-2, 0.1)
    p1 = heatmap(ε1, ε2, z; c=cgrad(:viridis, rev=true), clims, xlabel="ε1 and ε3", ylabel="ε2", title="log10(1-MPU)")
    foreach((c, levels) -> contour!(p1, ε1, ε2, contours; c, levels), cs, ls)

    z = map(log10 ∘ LD, sols |> permutedims)
    clims = (-2, 0.1)
    p2 = heatmap(ε1, ε2, z; c=cgrad(:viridis, rev=true), clims, xlabel="ε1 and ε3", ylabel="ε2", title="log10(LD)", cbar=false)
    foreach((c, levels) -> contour!(p2, ε1, ε2, contours; c, levels), cs, ls)

    z = map(x -> x.excgap, sols |> permutedims)
    clims = (0, 0.3)
    p3 = heatmap(ε1, ε2, z; c=cgrad(:viridis, rev=false), clims, xlabel="ε1 and ε3", ylabel="ε2", title="excgap")
    foreach((c, levels) -> contour!(p2, ε1, ε2, contours; c, levels), cs, ls)
    return plot(p2, p3; layout=(1, 2), size=0.7 .* (1200, 400))
end
mapslices(display ∘ slice_plot, sols3d, dims=[2, 3])
##

map(MPU, sols3d[1, :, :])

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
plot(map(x -> x.excgap, sols2))
plot!(map(x -> x.excgap, sols3))
plot(map(x -> x.gap, sols2))
plot!(map(x -> x.gap, sols3))

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