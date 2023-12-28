using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra#, BlackBoxOptim
using Plots
using Symbolics
using Folds
using Accessors
# using Optimization, OptimizationBBO, OptimizationNLopt, OptimizationOptimJL, OptimizationNOMAD, OptimizationEvolutionary, OptimizationMetaheuristics, OptimizationMultistartOptimization, OptimizationGCMAES, OptimizationCMAEvolutionStrategy, OptimizationPRIMA, Ipopt
function plotmajcoeffs(ws, zs)
    k = string.(keys(ws).values)
    labels = ["w", "z"] |> permutedims
    plot(k, [abs2.(ws).values, abs2.(zs).values], label=labels, legend=:topleft)
end

##
# N = 5
# c1 = FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
# c2 = FermionBdGBasis(1:N, (:↑, :↓))
# transport = Transport(QuantumDots.Pauli(), (; T=1 / 40, μ=(0.0, 0.0)))
##
N = 2
fixedparams = (; t=0.5, θ=parameter(π / 2, :diff), V=0, U=1 * 10, Ez=2)
c = FermionBasis(1:N, (:↑, :↓), qn=QuantumDots.parity)
cbdg = FermionBdGBasis(1:N, (:↑, :↓))
@variables Δ, εs[1:div(N + 1, 2)]::Real
@time symbolic_ham = hamiltonian(c; Δ, ε=reflect(εs, N), fixedparams...);
build_function(symbolic_ham, [Δ, εs...], expression=Val{false})
params = merge(fixedparams, (; Δ, ε=reflect(εs[1:1], N)))

f, f! = LongerPoorMansMajoranas.build_whamiltonian(c; params...)
cache = if c isa FermionBasis
    blockdiagonal(f([1.0im, 2.0]), c)
elseif c isa FermionBdGBasis
    f([1.0, 2.0]) |> Matrix |> BdGMatrix
end
@time f!(cache, [4.0, 2.1]);
##
@time f, f!, cache = hamfunc(RΔ_Rδϕ_Rε(), FermionBdGBasis(1:2, (:↑, :↓)), fixedparams);
h1 = f([1, 2, 0])
##
fp2 = @insert fixedparams.Δ = 1
@time f, f!, cache = hamfunc(Aϕ_Rε(), FermionBdGBasis(1:2, (:↑, :↓)), fp2);
h2 = f([0, 2, 0])
@assert norm(h1 - h2) < 1e-6
##
c = FermionBasis(1:3, (:↑, :↓), qn=QuantumDots.parity)
fixedparams = (; t=0.5, θ=parameter(2atan(2.0), :diff), V=0.2, U=2.5, Ez=1.25)
f, f!, cache = LongerPoorMansMajoranas.get_hamfuncs(c, fixedparams);
H = f([1, 1, 1])
tr(H)
eigvals(H)

sol = fullsolve(H, c)
##
c = FermionBasis(1:3, (:↑, :↓), qn=QuantumDots.parity)
fixedparams = (; t=0.5, θ=parameter(pi / 2, :diff), V=0, U=1, Ez=1.5)
prob = OptProb(; basis=c, fixedparams)
sols = [solve(prob, alg; maxiters=10000, MaxTime=10) for alg in best_algs()]
@time solve(prob, best_algs()[1]; maxiters=10000, MaxTime=1);
@profview solve(prob, best_algs()[5]; maxiters=10000, MaxTime=5)
@time solve(prob, best_algs()[5]; maxiters=10000, MaxTime=5);
bar(best_alg_names(), map(x -> LDf(x.optsol), sols); xrotation=30)
bar(best_alg_names(), map(x -> MPU(x.optsol), sols); xrotation=30)
bar(best_alg_names(), map(x -> MPU(x.optsol) + LDf(x.optsol), sols); xrotation=30)
##
Us = collect(range(0.0, 5, length=6))
Ezs = collect(range(0.5, 2.5, length=7))
itr = Base.product(Us, Ezs) |> collect
fixedparams = (; t=0.5, θ=parameter(pi / 2, :diff), V=0, Δ=1)
Nsols = Vector{Any}(undef, 3)
Threads.@threads for (n, N) in collect(enumerate(2:4))
    c = FermionBasis(1:N, (:↑, :↓), qn=QuantumDots.parity)
    UEsols = Matrix{Any}(undef, length(Us), length(Ezs))
    for (i, U) in enumerate(Us), (k, Ez) in enumerate(Ezs)
        prob = OptProb(; basis=c, fixedparams=merge(fixedparams, (; U, Ez)))
        algsols = [solve(prob, alg; maxiters=1000, MaxTime=5) for alg in best_algs()]
        UEsols[i, k] = algsols
    end
    Nsols[n] = UEsols
end
##
heatmap(Us, Ezs, map(x -> minimum(map(y -> MPU(y.optsol), x)), Nsols[1])')
heatmap(Us, Ezs, map(x -> minimum(map(y -> MPU(y.optsol), x)), Nsols[2])')
heatmap(Us, Ezs, map(x -> minimum(map(y -> MPU(y.optsol), x)), Nsols[3])')
map(x -> findmin(map(y -> MPU(y.optsol), x)), Nsols[1])
map(x -> findmin(map(y -> MPU(y.optsol), x)), Nsols[2])
map(x -> findmin(map(y -> MPU(y.optsol), x)), Nsols[3])
best_algs()[5]
##
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1)
hamfunc = LongerPoorMansMajoranas.hamfunc_ϕε
basis = FermionBasis(1:2, (:↑, :↓), qn=QuantumDots.parity)
basis = FermionBdGBasis(1:2, (:↑, :↓))
# algs = best_algs()[[1, 5]]
# algs = best_algs()[[1]]
algs = best_algs()[[1]]
UEsols2 = let Us = collect(range(0.0, 6, length=5)), Ezs = collect(range(0.1, 10, length=8))
    UEsols = Matrix{Any}(undef, length(Us), length(Ezs))
    for (i, U) in collect(enumerate(Us))
        for (k, Ez) in enumerate(Ezs)
            prob = OptProb(; hamfunc, basis, fixedparams=merge(fixedparams, (; U, Ez)))
            algsols = [solve(prob, alg; minexcgap=0.2, maxiters=1000, MaxTime=1) for alg in algs]
            UEsols[i, k] = algsols
        end
    end
    UEsols
end
##
fixedparams = (; t=0.5, θ=parameter(2atan(1.1), :diff), V=0, Δ=0.1, U=0, Ez=2.5)
smb = let basis = FermionBasis(1:2, (:↑, :↓), qn=QuantumDots.parity), ps = [1, 2]
    f, f!, c = LongerPoorMansMajoranas.hamfunc_ϕε(basis, fixedparams)
    f!(c, ps)
    fullsolve(c, basis)
end;
sbdg = let basis = FermionBdGBasis(1:2, (:↑, :↓)), ps = [1, 2]
    f, f!, c = LongerPoorMansMajoranas.hamfunc_ϕε(basis, fixedparams)
    h = f(ps)
    h = f!(c, ps)
    # f(ps) - c
    fullsolve(c, basis)
end;
smb.excgap
sbdg.excgap
LDf.([smb, sbdg])
##
# basis = FermionBasis(1:3, (:↑, :↓), qn=QuantumDots.parity)
basis = FermionBdGBasis(1:2, (:↑, :↓))
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1.0, U=0, Ez=4.0)
# optparams = Aϕ_Rε()
optparams = Rδϕ_Rε()
f, f!, cache = hamfunc(optparams, basis, fixedparams)
prob1 = OptProb(; hamfunc=x -> f!(cache, x), basis, optparams)
sol1 = solve(prob1, best_algs()[1]; minexcgap=0.2, maxiters=100000, MaxTime=10)
sol1.optsol |> LD
sol1.optsol |> MPU
##
function longersol(sol, fixedparams, N)
    δϕ = sol[1]
    ϕ = δϕ .* (0:N-1)
    ε = parameter(sol[2], :homogeneous)
    Δ = @. exp(1im * ϕ) * fixedparams.Δ
    newparams = merge(fixedparams, (; Δ, ε))
    basis = FermionBdGBasis(1:N, (:↑, :↓))
    h = LongerPoorMansMajoranas.whamiltonian(basis; newparams...) |> hermitianpart! |> BdGMatrix
    fullsolve(h, basis)
end
sols = [longersol(sol1.sol, n) for n in 2:10];
MPU.(sols)
##
f, f!, cache = hamfunc(Aϕ_Rε(), basis, fixedparams)

prob2 = OptProb(; hamfunc=LongerPoorMansMajoranas.hamfunc_rϕε, basis, fixedparams=merge(fixedparams, (; U=1, Ez=2)))
sol2 = solve(prob2, best_algs()[1]; minexcgap=0.3, maxiters=10000, MaxTime=5)
##
sol = solve(prob, best_algs()[7]; minexcgap=0.4, maxiters=1000, MaxTime=1)
sol.optsol.excgap
sol.optsol.gap
##
heatmap(map(x -> minimum(map(y -> log(MPU(y.optsol)), x)), UEsols2)'; c=cgrad(:viridis, rev=true), clims=(-6, 0))
heatmap(map(x -> minimum(map(y -> log(LD(y.optsol)), x)), UEsols2)'; c=cgrad(:viridis, rev=true), clims=(-6, 0))
heatmap(map(x -> x[1].optsol.excgap, UEsols)', clims=(0, 0.5))
heatmap(map(x -> tanh(x[1].optsol.gap), UEsols)', clims=(-0.01, 0.01), c=:redsblues)

heatmap(map(x -> minimum(map(y -> LD(y.optsol), x)), UEsols)')
##
map(x -> findmin(y -> log(LD(y.optsol)), x)[2], UEsols)
##
map(x -> x.params, sols)


##
Δ0 = 10
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=Δ0, U=0.0)
N = 2
basis1 = FermionBdGBasis(1:N, (:↑, :↓))
basis2 = basis1 #FermionBasis(1:N, (:↑, :↓), qn=QuantumDots.parity)
op1 = Rδϕ_Rε()
op2 = Aϕ_Rε()
alginds = [1, 1, 2, 3, 4]
Ezs_big = collect(range(0.1, 8, length=50)) .* Δ0
Ezsols1_big, Ezsols2 = let Ezs = Ezs_big
    Ezsols1 = Vector{Any}(undef, length(Ezs))
    Ezsols2 = Vector{Any}(undef, length(Ezs))
    for (k, Ez) in collect(enumerate(Ezs))
        fp = merge(fixedparams, (; Ez))
        f1, f1!, cache1 = hamfunc(op1, basis1, fp)
        # f2, f2!, cache2 = hamfunc(op2, basis2, fp)
        prob1 = OptProb(; hamfunc=x -> f1!(cache1, x), basis=basis1, optparams=op1)
        # prob2 = OptProb(; hamfunc=f2, basis=basis2, optparams=op2)
        algsols1 = [solve(prob1, alg; minexcgap=0.2, maxiters=10000, MaxTime=5) for alg in best_algs()[alginds]]
        # algsols2 = [solve(prob2, alg; minexcgap=0.2, maxiters=10000, MaxTime=2) for alg in best_algs()[alginds]]
        Ezsols1[k] = algsols1
        # Ezsols2[k] = algsols2
    end
    Ezsols1, Ezsols2
end
##
ldata = wload(datadir("Ezsols_big.jld2")) #todo for tomorrow
wsave(datadir("Ezsols_big.jld2"), Dict(zip(string.(eachindex(Ezsols1_big)), Ezsols1_big)))
##
bestsols = map(x -> findmin(y -> MPU((y.optsol)), x)[2], Ezsols1)
bestsols_big = map(x -> findmin(y -> MPU((y.optsol)), x)[2], Ezsols1_big)
result = Dict("optparams"=>map((s,n) -> collect(s[n].sol), Ezsols1_big, bestsols_big), "Ez" => Ezs_big, "fixedparams" => fixedparams)
wsave(datadir("Ezsols_big.jld2"), result)
##
Ns = 15:25
sols = [longersol(Ezsols1[k][bestsols[k]].sol, merge(fixedparams, (; Ez)), N) for N in Ns, (k, Ez) in enumerate(Ezs)];
sols_big = [longersol(Ezsols1_big[k][bestsols_big[k]].sol, merge(fixedparams, (; Ez)), N) for N in Ns, (k, Ez) in enumerate(Ezs_big)];
sols_big = [longersol(ldata["optparams"][k], merge(fixedparams, (; Ez)), N) for N in Ns, (k, Ez) in enumerate(ldata["Ez"])];
Ezs_big = ldata["Ez"]
##
heatmap(Ezs, Ns, map(log ∘ MPU, sols), c=cgrad(:viridis, rev=true), clims=(-6, 0), xlabel="Ez/Δ", ylabel="N", title="log(1-MPU)")
pmpu = heatmap(Ezs_big ./ Δ0, Ns, map(log ∘ MPU, sols_big), c=cgrad(:viridis, rev=true), clims=(-6, 0), xlabel="Ez/Δ", ylabel="N", title="log(1-MPU)")
##
heatmap(Ezs, Ns, map(log ∘ LD, sols), c=cgrad(:viridis, rev=true), clims=(-6, 0), xlabel="Ez/Δ", ylabel="N", title="log(LD)")
pld = heatmap(Ezs_big ./ Δ0, Ns, map(log ∘ LDf, sols_big), c=cgrad(:viridis, rev=true), clims=(-6, 0), xlabel="Ez/Δ", ylabel="N", title="log(LD)")
##
heatmap(Ezs, Ns, map(x -> log(abs(x.gap)), sols), c=cgrad(:viridis, rev=true), clims=(-5, 0), xlabel="Ez/Δ", ylabel="N", title="Egap")
pgap = heatmap(Ezs, Ns, map(x -> log(abs(x.gap)), sols_big), c=cgrad(:viridis, rev=true), clims=(-5, 0), xlabel="Ez/Δ", ylabel="N", title="log(Egap)")
##
heatmap(Ezs, Ns, map(x -> log(abs(x.excgap)), sols), c=cgrad(:viridis, rev=true), clims=(-3, 0), xlabel="Ez/Δ", ylabel="N", title="Excgap")
pexcgap = heatmap(Ezs, Ns, map(x -> (abs(x.excgap)), sols_big), c=cgrad(:viridis, rev=true), clims=(0, 0.2), xlabel="Ez/Δ", ylabel="N", title="Excgap")
##
heatmap(Ezs, Ns, map(x -> log(abs(x.gapratio)), sols_big), c=cgrad(:viridis, rev=true), clims=(-20, 0), xlabel="Ez/Δ", ylabel="N", title="gapratio")
##
plot(pmpu, pld, pgap, pexcgap; plot_title="Δ = 20t", size=(800, 600))
##
[display(plotmajcoeffs(s.majcoeffs...)) for s in sols_big[[1, 5, 10, 39], 80]];
[display(plotmajcoeffs(s.majcoeffs...)) for s in sols_big[[1, 5, 10, 39], 52]];
## 
δϕ = π

plot(map(x -> minimum(map(y -> log(MPU(y.optsol)), x)), Ezsols1); clims=(-6, 0))
plot!(map(x -> minimum(map(y -> log(1e-10 + MPU(y.optsol)), x)), Ezsols2); clims=(-6, 0))
##
plot(map(x -> log(MPU(x[1].optsol)), Ezsols1); clims=(-6, 0));
plot!(map(x -> log(MPU(x[2].optsol)), Ezsols2); clims=(-6, 0));
plot!(map(x -> log(MPU(x[1].optsol)), Ezsols1); clims=(-6, 0));
plot!(map(x -> log(MPU(x[2].optsol)), Ezsols2); clims=(-6, 0))
##
plot(map(x -> log(abs(x[1].optsol.gap)), Ezsols1); clims=(-6, 0));
plot!(map(x -> log(abs(x[2].optsol.gap)), Ezsols2); clims=(-6, 0));
plot!(map(x -> log(abs(x[1].optsol.gap)), Ezsols1); clims=(-6, 0));
plot!(map(x -> log(abs(x[2].optsol.gap)), Ezsols2); clims=(-6, 0))
##
plot(map(x -> ((x[1].optsol.excgap)), Ezsols1); clims=(-6, 0));
plot!(map(x -> ((x[2].optsol.excgap)), Ezsols2); clims=(-6, 0));
plot!(map(x -> ((x[1].optsol.excgap)), Ezsols1); clims=(-6, 0));
plot!(map(x -> ((x[2].optsol.excgap)), Ezsols2); clims=(-6, 0))
##
plot(map(x -> minimum(map(y -> log(LD(y.optsol)), x)), Ezsols1); clims=(-6, 0));
plot!(map(x -> minimum(map(y -> log(LD(y.optsol)), x)), Ezsols2); clims=(-6, 0))
##
map(x -> x[2].optsol.gap, Ezsols1)
map(x -> x[2].optsol.excgap, Ezsols1)

map(x -> x[2].sol |> collect, Ezsols1)
map(x -> x[2].sol |> collect, Ezsols2)
map(x -> x[2].sol[1:2] |> diff |> only, Ezsols2)
bestols = map(x -> findmin(y -> MPU((y.optsol)), x)[2], Ezsols1)


##
multisolve(f, alg, basis; MaxTime=5, minexcgap=1 / 4, exps=collect(range(0.1, 3, length=4)), maxiters=1000, initials=get_initials(basis), kwargs...) = multisolve(f, (alg,), basis; MaxTime, minexcgap, exps, maxiters, initials, kwargs...)
multisolve_cons(f, alg, basis; MaxTime=5, minexcgap=1 / 4, exps=collect(range(0.1, 3, length=4)), maxiters=1000, initials=get_initials(basis), kwargs...) = multisolve2(f, (alg,), basis; MaxTime, minexcgap, exps, maxiters, initials, kwargs...)

##
optimalgs = [ParticleSwarm(), SAMIN()]
nomadalgs = [NOMADOpt()]
bboalgs = [BBO_dxnes(), BBO_xnes(), BBO_de_rand_1_bin_radiuslimited(), BBO_adaptive_de_rand_1_bin_radiuslimited(), BBO_generating_set_search(), BBO_adaptive_de_rand_1_bin(), BBO_separable_nes(), BBO_probabilistic_descent(), BBO_resampling_memetic_search(), BBO_resampling_inheritance_memetic_search()]
bboalgs2 = [BBO_dxnes(), BBO_xnes(), BBO_de_rand_1_bin_radiuslimited(), BBO_adaptive_de_rand_1_bin_radiuslimited(), BBO_generating_set_search(), BBO_adaptive_de_rand_1_bin(), BBO_separable_nes(), BBO_probabilistic_descent()]
evoalgs = [OptimizationEvolutionary.GA(), OptimizationEvolutionary.DE(), ES(), CMAES()]
nloalgs = [NLopt.GN_DIRECT(), NLopt.GN_DIRECT_L(), NLopt.GN_DIRECT_L_RAND(), NLopt.GN_DIRECT_NOSCAL(), NLopt.GN_DIRECT_L_NOSCAL(), NLopt.GN_DIRECT_L_RAND_NOSCAL(), NLopt.GD_STOGO(), NLopt.GD_STOGO_RAND(), NLopt.GN_CRS2_LM(), NLopt.GN_ESCH()]
nlolocalgs = [NLopt.LN_PRAXIS(), NLopt.LN_COBYLA(), NLopt.LN_NEWUOA(), NLopt.LN_NEWUOA_BOUND(), NLopt.LN_NELDERMEAD(), NLopt.LN_SBPLX(), NLopt.LN_BOBYQA()]

msoalgs = [MultistartOptimization.TikTak(100)] # use in combination with local optimizer
metaalgs = [Metaheuristics.ECA(), Metaheuristics.DE(), PSO(), ABC(), CGSA(), SA(), WOA()]
metaalgs2 = [Metaheuristics.DE(), PSO(), CGSA(), SA(), WOA()]
gcmaesalgs = [GCMAESOpt()]
evo2algs = [CMAEvolutionStrategyOpt()]
primaalgs = [UOBYQA(), NEWUOA(), BOBYQA(), LINCOA(), COBYLA()]
primaalgs2 = [BOBYQA(), LINCOA()]
gradalgs = [Optim.NelderMead(), Optim.ConjugateGradient(), Optim.GradientDescent()]
consalgs = [Optim.IPNewton()]

##
c = FermionBasis(1:2, (:↑, :↓), qn=QuantumDots.parity)
fixedparams = (; t=0.5, θ=parameter(pi / 2, :diff), V=0, U=1, Ez=1.5)
f0, fs0 = opt_func(fixedparams, c; ad=nothing);
f, fs = opt_func2(fixedparams, c; ad=AutoFiniteDiff());
@time sol0 = multisolve(f0, BBO_probabilistic_descent(), c; maxiters=1000, MaxTime=2, exps=[0.1, 1]);
@time sol1 = multisolve2(f, Optim.IPNewton(), c; maxiters=1000, MaxTime=2, allow_f_increases=true, store_trace=false, g_tol=1e-3);
# @time sol2 = multisolve2(f, Optim.Newton(), c; maxiters=1000, MaxTime=10);
@time sol2 = multisolve2(f, Ipopt.Optimizer(), c; maxiters=1000, MaxTime=2);
# @time sol2 = multisolve2(f, Optim.IPNewton(), c; initials = sol.u, maxiters=100, MaxTime=10, exps=[1,2, 3]);
# @time sol2 = multisolve2(f, Optim.Newton(), c; initials = sol.u, maxiters=100, MaxTime=10, exps=[1,2,3]);
# @time sol2 = multisolve(f, (MultistartOptimization.TikTak(100), BBO_probabilistic_descent()), c, maxiters=10);

optsols = map(fs, [sol0, sol1, sol2]);
map(LDf, optsols)
map(x -> x.gap, optsols)
map(x -> x.excgap, optsols)
map(x -> x.excgap, optsols)
params = NamedTuple(zip((:Δ, :δϕ, :ε), decompose(sol)))
params2 = NamedTuple(zip((:Δ, :δϕ, :ε), decompose(sol2)))

##
results = []
c = FermionBasis(1:3, (:↑, :↓), qn=QuantumDots.parity)
fixedparams = (; t=0.5, θ=parameter(pi / 2, :diff), V=0, U=1, Ez=1.5)
f, fs = opt_func(fixedparams, c; ad=nothing);
fgrad, fsgrad = opt_func(fixedparams, c; ad=AutoFiniteDiff());
fgradcons, fsgradcons = opt_func_cons(fixedparams, c; ad=AutoFiniteDiff());
algs = vcat(bboalgs2, optimalgs, metaalgs2, evo2algs, nloalgs)
maxiters = 50
MaxTime = 5
for alg in algs
    println(alg)
    sol = multisolve(f0, alg, c; maxiters, MaxTime)
    optsol = fs(sol)
    params = NamedTuple(zip((:Δ, :δϕ, :ε), decompose(sol)))
    push!(results, (; alg, sol, optsol, params))
end
for alg in gradalgs
    println(alg)
    sol = multisolve(fgrad, alg, c; maxiters, MaxTime)
    optsol = fsgrad(sol)
    params = NamedTuple(zip((:Δ, :δϕ, :ε), decompose(sol)))
    push!(results, (; alg, sol, optsol, params))
end
for alg in consalgs
    println(alg)
    sol = multisolve_cons(fgradcons, alg, c; maxiters, MaxTime)
    optsol = fsgradcons(sol)
    params = NamedTuple(zip((:Δ, :δϕ, :ε), decompose(sol)))
    push!(results, (; alg, sol, optsol, params))
end
##
testedalgs = vcat(algs, gradalgs, consalgs)
bar(1:length(results), map(x -> x.optsol |> MPU, results); title="MPU")
bar(1:length(results), map(x -> x.optsol |> LDf, results); title="LDf")
bar(map(x -> x.sol.solve_time, results)[1:end-1]; title="time")
map(x -> x.sol[2], results) |> bar

# Bad times: 9,15.
# Kinda bad times: 12:16
##
push!(permsMPU, sortperm(map(x -> x.optsol |> MPU, results)))
push!(permsLDf, sortperm(map(x -> x.optsol |> LDf, results)))
##
permsMPU = []
permsLDf = []
##
optsols = [optimize(fixedparams, cbdg; MaxTime=2, tracemode=:compact) for N in 2:4];
# optsols2 = [optimize(N, fixedparams; MaxTime=N^2, initials=sol.bc, tracemode=:compact) for (N, sol) in zip(2:4, optsols)];
# optsols3 = [optimize(N, fixedparams; MaxTime=N^2, initials=sol.bc) for (N, sol) in zip(2:4, optsols)]
##
sols = optsols
for sol in sols
    plotmajcoeffs(sol.optsol.majcoeffs...) |> display
end
map(sol -> LDf(sol.optsol), sols)
map(sol -> MPU(sol.optsol), sols)
map(sol -> decompose(sol.bc), sols)
map(sol -> sol.optsol.excgap, sols)
##
for sol in optsols2
    plotmajcoeffs(sol.optsol.majcoeffs...) |> display
end
map(sol -> LDf(sol.optsol), optsols2)
map(sol -> MPU(sol.optsol), optsols2)
map(sol -> sol.bc, optsols2)
map(sol -> sol.optsol.excgap, optsols2)
## Reproduce old figure
Us = collect(range(0.0, 5, length=8))
Ezs = collect(range(0.5, 2.5, length=10))
itr = Base.product(Us, Ezs) |> collect
# fixedparams = (; t=0.5, θ=parameter(2atan(0.2), :diff), V=0)
fixedparams = (; t=0.5, θ=parameter(pi / 2, :diff), V=0)
c = FermionBasis(1:4, (:↑, :↓), qn=QuantumDots.parity)
sols = map((UEz) -> optimize(merge(fixedparams, (; U=UEz[1], Ez=UEz[2])), c; MaxTime=10, tracemode=:silent), itr);
##
heatmap(Us, Ezs, map(x -> MPU(x.optsol), sols)'; c=:viridis, clims=(0, 2))
heatmap(Us, Ezs, map(x -> x.optsol.excgap, sols)'; c=:viridis, clims=(0, 2))
heatmap(Us, Ezs, map(x -> x.params.Δ[2], sols)'; c=:viridis, clims=(0, 2))
# heatmap(last.(itr))
##
#check all methods in BlackBoxOptim.jl
fixedparams = (; t=0.5, θ=parameter(pi / 2, :diff), V=0)
c = FermionBasis(1:3, (:↑, :↓), qn=QuantumDots.parity)
# Methods = [:dxnes, :xnes, :de_rand_1_bin_radiuslimited, :adaptive_de_rand_1_bin_radiuslimited, :generating_set_search, :adaptive_de_rand_1_bin, :separable_nes, :probabilistic_descent, :resampling_memetic_search, :resampling_inheritance_memetic_search];
Methods = [:generating_set_search, :separable_nes, :probabilistic_descent, :resampling_inheritance_memetic_search];
Us = collect(range(0.5, 1.5, length=3))
Ezs = collect(range(0.9, 1.5, length=3))
optsols = Matrix{Any}(undef, length(Us), length(Ezs));
for (i, U) in enumerate(Us), (k, Ez) in enumerate(Ezs)
    optsols[i, k] = map(Method -> optimize(merge(fixedparams, (; U, Ez)), c; MaxTime=2, tracemode=:silent, Method), Methods)
end

##
optsols2 = Folds.map(Method -> optimize(merge(fixedparams, (; U=1.0, Ez=1.5)), c; MaxTime=1, tracemode=:silent, Method), Methods)
optsols3 = Folds.map(Method -> optimize(merge(fixedparams, (; U=1.0, Ez=1.5)), c; MaxTime=1, tracemode=:silent, Method), Methods)

##
optsols[findmin(x -> MPU(x.optsol), optsols)[2]].params
optsols2[findmin(x -> MPU(x.optsol), optsols2)[2]].params
optsols3[findmin(x -> MPU(x.optsol), optsols3)[2]].params
findmin(x -> MPU(x.optsol), optsols)
findmin(x -> MPU(x.optsol), optsols2)
findmin(x -> MPU(x.optsol), optsols3)

map(x -> MPU(x.optsol), optsols)

bar(string.(Methods), map(x -> MPU(x.optsol), optsols3), size=(1500, 500))