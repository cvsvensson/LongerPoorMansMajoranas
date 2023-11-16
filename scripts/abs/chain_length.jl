using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra#, BlackBoxOptim
using Plots
using Symbolics
using Folds
using Optimization, OptimizationBBO, OptimizationNLopt, OptimizationOptimJL, OptimizationNOMAD, OptimizationEvolutionary, OptimizationMetaheuristics, OptimizationMultistartOptimization, OptimizationGCMAES, OptimizationCMAEvolutionStrategy, OptimizationPRIMA
function plotmajcoeffs(ws, zs)
    k = string.(keys(ws).values)
    labels = ["w", "z"]
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
# L == 2N + div(N + 1, 2)
function decompose(ps, N=div(2length(ps), 3))#div(length(ps), 2))#div(length(ps),2)) 
    Nhalf = div(N + 1, 2)
    Nhalf2 = div(N, 2)
    rs = ps[1:Nhalf]
    ϕs = ps[Nhalf+1:Nhalf+Nhalf2]
    εs = ps[Nhalf+Nhalf2+1:end]
    return rs, ϕs, εs
end
function get_hamfuncs(c, fixedparams)
    N = div(QuantumDots.nbr_of_fermions(c), 2)
    Nhalf = div(N + 1, 2)
    @variables Δ[1:N], εs[1:Nhalf]::Real
    params = merge(fixedparams, (; Δ=collect(Δ[1:N]), ε=reflect(εs, N)))
    f, f! = LongerPoorMansMajoranas.build_whamiltonian(c; params...)
    fphase(Δs, εs) = f([Δs..., εs...])
    fphase(rs, δϕs, ϵs) = fphase((reflect(rs, N) .* exp.(1im .* diffreflect(δϕs, N))), ϵs)
    fphase(ps) = fphase(decompose(ps)...)
    fphase!(out, Δs, εs) = f!(out, [Δs..., εs...])
    fphase!(out, rs, δϕs, ϵs) = fphase!(out, (reflect(rs, N) .* exp.(1im .* diffreflect(δϕs, N))), ϵs)
    fphase!(out, ps) = fphase!(out, decompose(ps)...)
    ps = rand(2Nhalf + div(N, 2))
    cache = if c isa FermionBasis
        blockdiagonal(fphase(ps), c)
    elseif c isa FermionBdGBasis # Wait for update to QuantumDots
        fphase(ps) |> Matrix |> BdGMatrix
    end
    return fphase, fphase!, cache
end

@time f, f!, cache = get_hamfuncs(FermionBasis(1:2, (:↑, :↓), qn=QuantumDots.parity), fixedparams);
f(1:3)
##
function opt_func2(fixedparams, basis; ad=nothing)
    f0, f!, cache = get_hamfuncs(basis, fixedparams)
    target = x -> MPU(x) + LDf(x)
    extra_cost(a, b) = 0.0
    fullsolve2(x) = fullsolve(f!(cache, x), basis)
    # f2(x) = fullsolve(f0(x), basis)
    cost(x, p) = target(fullsolve2(x))
    function cons(res, x, minexcgap)
        sol = fullsolve2(x)
        res[1] = sol.gap
        res[2] = sol.excgap - minexcgap
    end
    if isnothing(ad)
        return OptimizationFunction(cost; cons), fullsolve2
    else
        return OptimizationFunction(cost, ad; cons), fullsolve2
    end
end
function opt_func(fixedparams, basis; ad=nothing)
    f0, f!, cache = get_hamfuncs(basis, fixedparams)
    target = x -> MPU(x) + LDf(x)
    extra_cost(a, b) = 0.0
    f2(x) = fullsolve(f!(cache, x), basis)
    # f2(x) = fullsolve(f0(x), basis)
    function f(x, (exp, minexcgap))
        sol = f2(x)
        LongerPoorMansMajoranas.cost_function(sol.gap, sol.excgap, target(sol); exp, minexcgap) + extra_cost(x, exp)
    end
    if isnothing(ad)
        return OptimizationFunction(f), f2
    else
        return OptimizationFunction(f, ad), f2
    end
end
multisolve(f, alg, basis; MaxTime=5, minexcgap=1 / 4, exps=collect(range(0.1, 3, length=4)), maxiters=1000, initials=[ones(div(div(length(basis), 2) + 1, 2))..., pi .* ones(div(div(length(basis), 2), 2))..., zeros(div(div(length(basis), 2) + 1, 2))...]) = multisolve(f, (alg,), basis; MaxTime, minexcgap, exps, maxiters, initials)
multisolve2(f, alg, basis; MaxTime=5, minexcgap=1 / 4, exps=collect(range(0.1, 3, length=4)), maxiters=1000, initials=[ones(div(div(length(basis), 2) + 1, 2))..., pi .* ones(div(div(length(basis), 2), 2))..., zeros(div(div(length(basis), 2) + 1, 2))...]) = multisolve2(f, (alg,), basis; MaxTime, minexcgap, exps, maxiters, initials)
function multisolve(f, alg::Tuple, basis; MaxTime=5, minexcgap=1 / 4, exps=collect(range(0.1, 3, length=4)), maxiters=1000, initials=[ones(div(div(length(basis), 2) + 1, 2))..., pi .* ones(div(div(length(basis), 2), 2))..., zeros(div(div(length(basis), 2) + 1, 2))...])
    N = div(length(basis), 2)
    Nhalf = div(N + 1, 2)
    refinements = length(exps)
    maxtime = MaxTime / refinements
    Δranges = [(0.1, 2.0) for i in 1:Nhalf]
    δϕranges = [(0.0, 2.0pi) for i in 1:div(N, 2)]
    εranges = [(-20.0, 20.0) for i in 1:Nhalf]
    ranges = [Δranges..., δϕranges..., εranges...]
    lb = map(first, ranges)
    ub = map(last, ranges)
    println("Initial point: ", initials)
    prob = OptimizationProblem(f, initials, (first(exps), minexcgap); lb=map(first, ranges), ub=map(last, ranges))
    sol = solve(prob, alg...; maxiters, maxtime)
    for (n, exp) in enumerate(Iterators.drop(exps, 1))
        newinitials = map(clamp, sol.u, lb, ub)
        println("$n, Sweet spot:", newinitials)
        prob = OptimizationProblem(f, initials, (exp, minexcgap); lb=map(first, ranges), ub=map(last, ranges))
        sol = solve(prob, alg...; maxiters, maxtime)
    end
    return sol
end
function multisolve2(f, alg, basis; MaxTime=5, minexcgap=1 / 4, maxiters=1000, initials=[ones(div(div(length(basis), 2) + 1, 2))..., pi .* ones(div(div(length(basis), 2), 2))..., zeros(div(div(length(basis), 2) + 1, 2))...], kwargs...)
    N = div(length(basis), 2)
    Nhalf = div(N + 1, 2)
    Δranges = [(0.1, 2.0) for i in 1:Nhalf]
    δϕranges = [(0.0, 2.0pi) for i in 1:div(N, 2)]
    εranges = [(-20.0, 20.0) for i in 1:Nhalf]
    ranges = [Δranges..., δϕranges..., εranges...]
    lb = map(first, ranges)
    ub = map(last, ranges)
    println("Initial point: ", initials)
    lcons = [0.0, 0.0]
    ucons = [0.0, Inf]
    prob = OptimizationProblem(f, initials, minexcgap; lb, ub, lcons, ucons, kwargs...)#, allow_f_increases = true, store_trace = true)
    sol = solve(prob, alg; maxiters, maxtime=MaxTime, abstol=1e-6, reltol=1e-6)
    return sol
end
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

##
c = FermionBasis(1:3, (:↑, :↓), qn=QuantumDots.parity)
fixedparams = (; t=0.5, θ=parameter(pi / 4, :diff), V=0, U=1, Ez=1.5)
f0, fs0 = opt_func(fixedparams, c; ad=nothing);
f, fs = opt_func2(fixedparams, c; ad=AutoFiniteDiff());
@time sol0 = multisolve(f0, BBO_probabilistic_descent(), c; maxiters=1000, MaxTime=10, exps=[0.1, 1]);
@time sol1 = multisolve2(f, Optim.IPNewton(), c; maxiters=1000, MaxTime=10, allow_f_increases=true, store_trace=false, g_tol = 1e-3);
# @time sol2 = multisolve2(f, Optim.Newton(), c; maxiters=1000, MaxTime=10);
@time sol2 = multisolve2(f, Ipopt.Optimizer(), c; maxiters=1000, MaxTime=10);
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
algs = vcat(bboalgs2, optimalgs, metaalgs2, evo2algs, primaalgs2, nloalgs, gradalgs)
for alg in algs#algs
    println(alg)
    sol = multisolve(f, alg, c; maxiters=20)
    optsol = fs(sol)
    params = NamedTuple(zip((:Δ, :δϕ, :ε), decompose(sol)))
    push!(results, (; alg, sol, optsol, params))
end
##
bar(1:length(results), map(x -> x.optsol |> MPU, results); title="MPU")
bar(1:length(results), map(x -> x.optsol |> LDf, results); title="LDf")
bar(map(x -> x.sol.solve_time, results); title="time")
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