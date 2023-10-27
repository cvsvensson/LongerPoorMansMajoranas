Base.@kwdef struct Optimizer{f,r,i,t,ec}
    hamfunc::f
    ranges::Vector{r}
    initials::Vector{i}
    MaxTime::Int = 10
    minexcgap::Float64 = 0.0
    exps::Vector{Float64} = Float64.(collect(range(0.5, 3; length=4)))
    target::t = LD
    tracemode::Symbol = :silent
    extra_cost::ec = (x...) -> 0
    Method::Symbol = :probabilistic_descent
    PopulationSize::Int = 100
    TargetFitness::Float64 = 0.0
end


cost_function(energies, reduced::Number; exp=12.0, minexcgap=0) = cost_reduced(reduced) + cost_energy(energies; exp, minexcgap)
cost_energy(energies; minexcgap=0, exp) = cost_gapratio(gapratio(energies...); exp) + ((excgap(energies...) - minexcgap) < 0 ? 1 + 10.0^exp*abs(excgap(energies...) - minexcgap) : 0)
cost_gapratio(gr; exp) = abs(gr) > 2 * 10.0^(-exp) ? 1.0 + 10^(exp) * abs2(gr) : abs2(gr)
cost_reduced(reduced) = reduced^2



tracemode(opt::Optimizer) = opt.tracemode
function cost(exp, opt::Optimizer)
    function _cost(args)
        sol = fullsolve(opt.hamfunc(args...))
        cost_function(sol.energies, opt.target(sol); exp, opt.minexcgap) + opt.extra_cost(args, exp)
    end
end

function get_sweet_spot(opt::Optimizer)
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
    opt = Optimizer(;
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