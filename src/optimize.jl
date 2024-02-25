abstract type AbstractOptParams end
struct Aϕ_Rε <: AbstractOptParams end
struct RΔ_Rδϕ_Rε <: AbstractOptParams end
struct Rδϕ_Rε <: AbstractOptParams end
struct Hδϕ_Hε <: AbstractOptParams end

Base.@kwdef struct BBOptimizer{f,r,i,t,ec,B}
    hamfunc::f
    ranges::Vector{r}
    initials::Vector{i}
    MaxTime::Float64 = 10.0
    minexcgap::Float64 = 0.0
    exps::Vector{Float64} = Float64.(collect(range(0.5, 3; length=4)))
    target::t = LD
    tracemode::Symbol = :silent
    extra_cost::ec = (x...) -> 0
    Method::Symbol = :probabilistic_descent
    PopulationSize::Int = 100
    TargetFitness::Float64 = 0.0
    basis::B
end


Base.@kwdef struct Optimizer{f,r,i,t,ec,B}
    hamfunc::f
    ranges::Vector{r}
    initials::Vector{i}
    MaxTime::Float64 = 10.0
    minexcgap::Float64 = 0.0
    exps::Vector{Float64} = Float64.(collect(range(0.5, 3; length=4)))
    target::t = LD
    tracemode::Symbol = :silent
    extra_cost::ec = (x...) -> 0
    Method::Symbol = :probabilistic_descent
    PopulationSize::Int = 100
    TargetFitness::Float64 = 0.0
    basis::B
end


cost_function(gap, excgap, reduced::Number; exp, minexcgap) = cost_reduced(reduced) + cost_energy(gap, excgap; exp, minexcgap)
cost_energy(gap, excgap; minexcgap, exp) = cost_gap(gap, exp) + cost_excgap(excgap, minexcgap, exp)
cost_excgap(excgap, minexcgap, exp) = ((excgap - minexcgap) < 0 ? 1.0 + 10.0^exp * abs(excgap - minexcgap) : 0.0)
# cost_gap(gap, exp) = abs(gap) > 2 * 10.0^(-exp) ? 1.0 + 10^(exp) * abs2(gap) : abs2(gap)
cost_gap(gap, exp) = 10.0^(exp) * abs2(gap)
cost_reduced(reduced) = reduced^2

best_algs() = [BBO_probabilistic_descent(), BBO_generating_set_search(), BBO_adaptive_de_rand_1_bin_radiuslimited(), Metaheuristics.SA(), Metaheuristics.DE(), Optim.IPNewton(), Optim.ConjugateGradient(), Optim.NelderMead()]
best_alg_names() = string.([:BBO_probabilistic_descent, :BBO_generating_set_search, :BBO_adaptive_de_rand_1_bin_radiuslimited, :Metaheuristics_SA, :Metaheuristics_DE, :Optim_IPNewton, :Optim_ConjugateGradient, :Optim_NelderMead])
function decompose_rϕε(ps, N=div(2length(ps), 3))
    Nhalf = div(N + 1, 2)
    Nhalf2 = div(N, 2)
    rs = ps[1:Nhalf]
    ϕs = ps[Nhalf+1:Nhalf+Nhalf2]
    εs = ps[Nhalf+Nhalf2+1:end]
    return rs, ϕs, εs
end
function decompose_ϕε(ps, N=length(ps))#div(length(ps), 2))
    # Nhalf = div(N + 1, 2)
    Nhalf2 = div(N, 2)
    ϕs = ps[1:Nhalf2]
    εs = ps[Nhalf2+1:end]
    return ϕs, εs
end
function decompose_allϕ_ε(ps, N)#div(length(ps), 2))
    ϕs = ps[1:N]
    εs = ps[N+1:end]
    return ϕs, εs
end
splatter_rϕε(ps, N) = splatter_rϕε(decompose_rϕε(ps)..., N)
function splatter_rϕε(rs, δϕs, ϵs, N)
    # N = div(2 * (length(rs) + length(δϕs) + length(ϵs)), 3)
    return vcat(reflect(rs, N) .* exp.(1im .* diffreflect(δϕs, N)), ϵs)
end

hamfunc(::RΔ_Rδϕ_Rε, c, fixedparams) = hamfunc_rϕε(c, fixedparams)
hamfunc(::Rδϕ_Rε, c, fixedparams) = hamfunc_ϕε(c, fixedparams)
hamfunc(::Aϕ_Rε, c, fixedparams) = hamfunc_allϕ_ε(c, fixedparams)
hamfunc(::Hδϕ_Hε, c, fixedparams) = hamfunc_Hδϕ_Hε(c, fixedparams)
function hamfunc_rϕε(c, fixedparams)
    N = div(QuantumDots.nbr_of_fermions(c), 2)
    Nhalf = div(N + 1, 2)
    @variables Δ[1:N], εs[1:Nhalf]::Real
    params = merge(fixedparams, (; Δ=collect(Δ[1:N]), ε=reflect(εs, N)))
    f, f! = LongerPoorMansMajoranas.build_whamiltonian(c; params...)
    f2(ps) = f(splatter_rϕε(ps, N))
    f2!(out, ps) = f!(out, splatter_rϕε(ps, N))
    ps = rand(2Nhalf + div(N, 2))
    cache = get_cache(c, f2(ps))
    return f2, f2!, cache
end
function hamfunc_ϕε(c, fixedparams)
    N = div(QuantumDots.nbr_of_fermions(c), 2)
    Nhalf = div(N + 1, 2)
    Nhalf2 = div(N, 2)
    @variables Δ[1:N], εs[1:Nhalf]::Real
    params = merge(fixedparams, (; Δ=collect(Δ[1:N]), ε=reflect(εs, N)))
    f, f! = LongerPoorMansMajoranas.build_whamiltonian(c; params...)
    splatter_ϕε(ps, N) = splatter_ϕε(decompose_ϕε(ps, N)..., N)
    splatter_ϕε(δϕs, ϵs, N) = vcat(fixedparams.Δ .* exp.(1im .* diffreflect(δϕs, N)), ϵs)
    f2(ps) = f(splatter_ϕε(ps, N))
    f2!(out, ps) = f!(out, splatter_ϕε(ps, N))
    ps = rand(N)
    cache = get_cache(c, Matrix(f2(ps)))
    # display(f2(ps))
    return f2, f2!, cache
end
function hamfunc_Hδϕ_Hε(c, fixedparams)
    N = div(QuantumDots.nbr_of_fermions(c), 2)
    @variables Δ[1:N], ε::Real
    params = merge(fixedparams, (; Δ=collect(Δ[1:N]), ε))
    f, f! = LongerPoorMansMajoranas.build_whamiltonian(c; params...)
    function splatter(δϕϵs, N)
        ϕs = (0:N-1) .* (δϕϵs[1])
        vcat(fixedparams.Δ .* exp.(1im .* ϕs), δϕϵs[2])
    end
    f2(ps) = f(splatter(ps, N))
    f2!(out, ps) = f!(out, splatter(ps, N))
    ps = rand(2)
    cache = get_cache(c, Matrix(f2(ps)))
    # display(f2(ps))
    return f2, f2!, cache
end
function hamfunc_allϕ_ε(c, fixedparams)
    N = div(QuantumDots.nbr_of_fermions(c), 2)
    Nhalf = div(N + 1, 2)
    # Nhalf2 = div(N, 2)
    @variables Δ[1:N], εs[1:Nhalf]::Real
    params = merge(fixedparams, (; Δ=collect(Δ[1:N]), ε=reflect(εs, N)))
    f, f! = LongerPoorMansMajoranas.build_whamiltonian(c; params...)
    splatter_allϕ_ε(ps, N) = splatter_allϕ_ε(decompose_allϕ_ε(ps, N)..., N)
    splatter_allϕ_ε(ϕs, ϵs, N) = vcat(fixedparams.Δ .* exp.(1im .* ϕs), ϵs)
    f2(ps) = f(splatter_allϕ_ε(ps, N))
    f2!(out, ps) = f!(out, splatter_allϕ_ε(ps, N))
    ps = rand(N + Nhalf)
    cache = get_cache(c, f2(ps))
    return f2, f2!, cache
end



get_cache(c::FermionBasis, ham) = blockdiagonal(ham, c)
get_cache(c::FermionBdGBasis, ham) = (Matrix(ham))
# function get_hamfuncs(c, fixedparams, hamfunc)
# N = div(QuantumDots.nbr_of_fermions(c), 2)
# Nhalf = div(N + 1, 2)
# @variables Δ[1:N], εs[1:Nhalf]::Real
# params = merge(fixedparams, (; Δ=collect(Δ[1:N]), ε=reflect(εs, N)))
# f, f! = LongerPoorMansMajoranas.build_whamiltonian(c; params...)
# f2(ps) = f(splatter(ps))
# f2!(out, ps) = f!(out, splatter(ps))
# fphase(Δs, εs) = f([Δs..., εs...])
# fphase(rs, δϕs, ϵs) = fphase((reflect(rs, N) .* exp.(1im .* diffreflect(δϕs, N))), ϵs)
# fphase(ps) = fphase(decompose(ps)...)
# fphase!(out, Δs, εs) = f!(out, [Δs..., εs...])
# fphase!(out, rs, δϕs, ϵs) = fphase!(out, (reflect(rs, N) .* exp.(1im .* diffreflect(δϕs, N))), ϵs)
# fphase!(out, ps) = fphase!(out, decompose(ps)...)
# ps = rand(2Nhalf + div(N, 2))
#     cache = if c isa FermionBasis
#         blockdiagonal(f2(ps), c)
#     elseif c isa FermionBdGBasis # Wait for update to QuantumDots
#         f2(ps) |> Matrix |> BdGMatrix
#     end
#     return fphase, fphase!, cache
# end

@kwdef struct OptProb{F,B,P,T}
    hamfunc::F
    basis::B
    optparams::P
    target::T = x -> MPU(x) + LDf(x)
end

opt_func(opt::OptProb, ::IPNewton) = opt_func_cons(opt.hamfunc, opt.basis, AutoFiniteDiff(), opt.target)
opt_func(opt::OptProb, ::Optim.ConjugateGradient) = opt_func(opt.hamfunc, opt.basis, AutoFiniteDiff(), opt.target)
opt_func(opt::OptProb, ::NelderMead) = opt_func(opt.hamfunc, opt.basis, AutoFiniteDiff(), opt.target)
opt_func(opt::OptProb, alg) = opt_func(opt.hamfunc, opt.basis, nothing, opt.target)


function opt_func_cons(hamfunc, basis, ad=nothing, target=x -> MPU(x) + LDf(x))
    # f0, f!, cache = hamfunc(basis, fixedparams)
    extra_cost(a, b) = 0.0
    fullsolve2(x) = fullsolve(hamfunc(x), basis)
    cost(x, p) = target(fullsolve2(x))
    function cons(res, x, minexcgap)
        sol = fullsolve2(x)
        res[1] = sol.gap
        res[2] = min(sol.excgap - minexcgap, 0.0)
    end
    if isnothing(ad)
        return OptimizationFunction(cost; cons), fullsolve2
    else
        return OptimizationFunction(cost, ad; cons), fullsolve2
    end
end
function opt_func(hamfunc, basis, ad=nothing, target=x -> MPU(x) + LDf(x))
    # f0, f!, cache = hamfunc(basis, fixedparams)
    # target = x -> MPU(x) + LDf(x)
    extra_cost(a, b) = 0.0
    f2(x) = fullsolve(hamfunc(x), basis)
    # f2(x) = fullsolve(f0(x), basis)
    function f(x, (exp, minexcgap))
        sol = f2(x)
        cost_function(sol.gap, sol.excgap, target(sol); exp, minexcgap) + extra_cost(x, exp)
    end
    if isnothing(ad)
        return OptimizationFunction(f), f2
    else
        return OptimizationFunction(f, ad), f2
    end
end

get_initials(prob::OptProb) = get_initials(prob.optparams, div(length(prob.basis), 2))
get_ranges(prob::OptProb) = get_ranges(prob.optparams, div(length(prob.basis), 2))
initial_Aϕ(N) = zeros(N)
initial_Rε(N) = zeros(div(N + 1, 2))
initial_Hε(N) = [0.0]
initial_Hδϕ(N) = [0.0]
initial_Rδϕ(N) = 0.0 .* ones(div(N, 2))
initial_RΔ(N) = ones(div(N + 1, 2))
ranges_Aϕ(N) = [(0.0, 2.0pi) for i in 1:N]
ranges_Rε(N) = [100 .* (0, 1) for i in 1:div(N + 1, 2)]
ranges_Rδϕ(N) = [(0.0, 1.0pi) for i in 1:div(N, 2)]
ranges_RΔ(N) = [(0.01, 10.0) for i in 1:div(N + 1, 2)]
ranges_Hε(N) = [100 .* (0, 1)]
ranges_Hδϕ(N) = [(0.0, 1.0pi)]

function get_initials(::Aϕ_Rε, N)
    vcat(initial_Aϕ(N), initial_Rε(N))
end
function get_initials(::RΔ_Rδϕ_Rε, N)
    vcat(initial_RΔ(N), initial_Rδϕ(N), initial_Rε(N))
end
function get_initials(::Hδϕ_Hε, N)
    vcat(initial_Hδϕ(N), initial_Hε(N))
end
function get_initials(::Rδϕ_Rε, N)
    vcat(initial_Rδϕ(N), initial_Rε(N))
end
function get_ranges(::Rδϕ_Rε, N)
    δϕranges = ranges_Rδϕ(N)
    εranges = ranges_Rε(N)
    vcat(δϕranges, εranges)
end
function get_ranges(::Hδϕ_Hε, N)
    δϕranges = ranges_Hδϕ(N)
    εranges = ranges_Hε(N)
    vcat(δϕranges, εranges)
end
function get_ranges(::RΔ_Rδϕ_Rε, N)
    Δranges = ranges_RΔ(N)
    δϕranges = ranges_Rδϕ(N)
    εranges = ranges_Rε(N)
    vcat(Δranges, δϕranges, εranges)
end
function get_ranges(::Aϕ_Rε, N)
    ϕranges = ranges_Aϕ(N)
    εranges = ranges_Rε(N)
    vcat(ϕranges, εranges)
end
default_exps() = collect(range(0.1, 3, length=5))
function SciMLBase.solve(prob::OptProb, alg; MaxTime=5, minexcgap=1 / 4, exps=default_exps(), maxiters=1000, initials=get_initials(prob), kwargs...)
    f, fs = opt_func(prob, alg)
    refinements = length(exps)
    maxtime = MaxTime / refinements
    ranges = get_ranges(prob)
    lb = map(first, ranges)
    ub = map(last, ranges)
    newinitials = map(clamp, initials, lb, ub)
    println("Initial point: ", newinitials)
    prob2 = OptimizationProblem(f, newinitials, (first(exps), minexcgap); lb, ub)
    sol = solve(prob2, alg; maxiters, maxtime)
    for (n, exp) in enumerate(Iterators.drop(exps, 1))
        newinitials = map(clamp, sol.u, lb, ub)
        println("$n, Sweet spot:", newinitials)
        prob2 = OptimizationProblem(f, initials, (exp, minexcgap); lb=map(first, ranges), ub=map(last, ranges))
        sol = solve(prob2, alg; maxiters, maxtime, kwargs...)
    end
    optsol = fs(sol)
    # params = NamedTuple(zip((:Δ, :δϕ, :ε), decompose(sol)))
    return (; alg, sol, optsol)
end

function SciMLBase.solve(prob::OptProb, alg::IPNewton; MaxTime=5, minexcgap=1 / 4, maxiters=1000, initials=get_initials(prob), kwargs...)
    ranges = get_ranges(prob)
    lb = map(first, ranges)
    ub = map(last, ranges)
    println("Initial point: ", initials)
    lcons = [0.0, 0.0]
    ucons = [0.0, Inf]
    f, fs = opt_func(prob, alg)
    prob2 = OptimizationProblem(f, initials, minexcgap; lb, ub, lcons, ucons, kwargs...)#, allow_f_increases = true, store_trace = true)
    sol = solve(prob2, alg; maxiters, maxtime=MaxTime, kwargs...)
    optsol = fs(sol)
    # params = NamedTuple(zip((:Δ, :δϕ, :ε), decompose(sol)))
    return (; alg, sol, optsol)
    # return sol
end

tracemode(opt::BBOptimizer) = opt.tracemode
function cost(exp, opt::BBOptimizer)
    function _cost(args)
        sol = fullsolve(opt.hamfunc(args...), opt.basis)
        cost_function(sol.gap, sol.excgap, opt.target(sol); exp, minexcgap=opt.minexcgap) + opt.extra_cost(args, exp)
    end
end

function get_sweet_spot(opt::BBOptimizer)
    refinements = length(opt.exps)
    SearchRange = opt.ranges #map(expand_searchrange, opt.ranges, opt.initials)
    NumDimensions = length(SearchRange)
    MaxTime = opt.MaxTime / refinements
    TargetFitness = opt.TargetFitness
    println("Initial point: ", opt.initials)
    println("SearchRange: ", SearchRange)
    Method = opt.Method
    res = bboptimize(cost(first(opt.exps), opt), opt.initials;
        SearchRange, NumDimensions, MaxTime, TraceInterval=10.0,
        TraceMode=tracemode(opt),
        Method, PopulationSize=opt.PopulationSize,
        TargetFitness)
    for exp in Iterators.drop(opt.exps, 1)
        ss::typeof(opt.initials) = best_candidate(res)
        println("Sweet spot:", ss)
        println("SearchRange:", SearchRange)
        res = bboptimize(cost(exp, opt), ss; SearchRange, NumDimensions,
            MaxTime, TraceInterval=10.0, TraceMode=tracemode(opt),
            Method, TargetFitness,
            PopulationSize=opt.PopulationSize)
    end
    return res
end

##
function anti_parallel_sweet_spot(; Δ, tratio, h, U, V, t, MaxTime, exps=collect(range(0.5, 3, length=4)), target, kwargs...)
    fixedparams = (; Δ, tratio, h, U, V, t)
    μ1, μ2 = kitaev_μ_zero(Δ, h, U)
    ϕ = 0.5 * pi
    hamfunc(ϕ, μ1, μ2) = hamiltonian(c; μ1, μ2, ϕ, fixedparams...)
    maxh = max(abs.(h)...)
    maxU = max(abs.(U)...)
    opt = BBOptimizer(;
        hamfunc,
        ranges=[(0.0, 1.0π), (0.0, 1.1 * μ1 + maxh + maxU), (-maxh - maxU + μ2, maxU + V)],
        initials=Float64.([ϕ, μ1, μ2]),
        MaxTime, exps, target,
        tracemode=:silent,
        extra_cost=((ϕ, μ1, μ2), e) -> exp(-(e * abs(μ1 - μ2) + 1)^4) + e * (μ2 > μ1) * (μ2 - μ1),
        kwargs...)
    ss = best_candidate(get_sweet_spot(opt))
    optsol = solve(opt.hamfunc(ss...))
    parameters = merge(fixedparams, NamedTuple(zip((:ϕ, :μ1, :μ2), ss)))
    optimization = opt
    sweet_spot = merge(optsol, (; parameters, optimization))
    return sweet_spot
end

function sweet_spot_scan((xs, xlabel), (ys, ylabel), get_ss=anti_parallel_sweet_spot; fixedparams, MaxTime, target, Method=:adaptive_de_rand_1_bin_radiuslimited, kwargs...)
    iter = collect(Base.product(xs, ys))
    ss = Folds.map((xy) -> get_ss(; fixedparams..., Dict(xlabel => xy[1], ylabel => xy[2])..., MaxTime, target, Method, kwargs...), iter)
    return Dict(:sweet_spots => ss, :x => xs, :y => ys, :xlabel => xlabel, :ylabel => ylabel,
        :fixedparams => fixedparams, :MaxTime => MaxTime, :target => target, :Method => Method)
end