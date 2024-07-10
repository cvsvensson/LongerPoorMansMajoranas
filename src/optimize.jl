
@kwdef struct OptProb{F,B,P,T}
    hamfunc::F
    basis::B
    optparams::P
    target::T
end

struct BestOf{OS}
    optimizers::OS
    function BestOf(opts::OS) where {OS}
        if any(isa.(opts, Optim.ConjugateGradient))
            throw(ArgumentError("ConjugateGradient is not supported in BestOf"))
        end
        new{OS}(opts)
    end
end
struct NoPenalty end
(::NoPenalty)(x...) = 0
# struct OptProb{FS,T,PF}
#     fullsolve::FS
#     target::T
#     penalty::PF
#     function OptProb(fullsolve::FS, target::T, penalty::PF=NoPenalty()) where {FS,T,PF}
#         new{FS,T,PF}(fullsolve, target, penalty)
#     end
# end
struct NLOptProb{FS,T,C,P}
    eigfunc::FS
    target::T
    constraint::C
    params::P
end

struct ScheduledOptProb{FS,T,PF}
    eigfunc::FS
    target::T
    penalty::PF
    function ScheduledOptProb(eigfunc::FS, target::T, penalty::PF=NoPenalty()) where {FS,T,PF}
        new{FS,T,PF}(eigfunc, target, penalty)
    end
end
get_initials(prob::ScheduledOptProb) = get_initials(prob.prob)
get_ranges(prob::ScheduledOptProb) = get_ranges(prob.prob)
abstract type AbstractPenalty end
struct ScheduledPenalty{EC} <: AbstractPenalty
    extra_cost::EC # (sol, x, i) -> cost
end
struct ConstantPenalty{EC} <: AbstractPenalty
    extra_cost::EC # (sol, x) -> cost
end
(pf::ConstantPenalty)(sol, x, i) = pf.extra_cost(sol, x)
(pf::ScheduledPenalty)(sol, x, i) = pf.extra_cost(sol, x, i)
function GapPenalty(exps)
    ScheduledPenalty((sol, x, i) -> cost_gap(get_gap(sol), exps[i]))
end
function MinExcGapPenalty(minexcgap, exps)
    ScheduledPenalty((sol, x, i) -> cost_excgap(sol.excgap, minexcgap, exps[i]))
end
function Base.:+(p1::ScheduledPenalty, p2::ScheduledPenalty)
    ScheduledPenalty((sol, x, i) -> p1.extra_cost(sol, x, i) + p2.extra_cost(sol, x, i))
end

function SciMLBase.solve(prob::ScheduledOptProb, alg; kwargs...)
    function f(x, (i,))
        sol = prob.eigfunc(x)
        prob.target(sol) +
        prob.penalty(sol, x, i)
    end
    of = OptimizationFunction(f, AutoFiniteDiff())
    sol = __solve(of, alg; kwargs...)
    optsol = all_info(prob.eigfunc(sol.u))
    gap_der = get_gap_derivatives(prob.eigfunc, sol.u)
    return (; sol, optsol, gap_der...)
end

function SciMLBase.solve(prob::ScheduledOptProb, alg::BestOf; MaxTime, kwargs...)
    function f(x, (i,))
        sol = prob.eigfunc(x)
        prob.target(sol) +
        prob.penalty(sol, x, i)
    end
    of = OptimizationFunction(f, AutoFiniteDiff())
    MaxTime = MaxTime / length(alg.optimizers)
    _sols = map(alg -> __solve(of, alg; MaxTime, kwargs...), alg.optimizers)
    _res = map(sol -> (; sol, optsol=all_info(prob.eigfunc(sol)), get_gap_derivatives(prob.eigfunc, sol)...), _sols)
    # return _sols
    sols = sort(_res, by=x -> prob.target(x.optsol))
    return merge(sols[1], (; all_sols=sols))
end

function _default_nl_opt(prob::NLOptProb)
    alg = get(prob.params, :alg, :LN_COBYLA)
    opt = Opt(alg, prob.params.length)
    function f(x, grad)
        sol = prob.eigfunc(x)
        prob.target(sol)
    end
    function g(x, grad)
        sol = prob.eigfunc(x)
        prob.constraint(sol, x)
    end
    opt.min_objective = f
    opt.ftol_rel = get(prob.params, :ftol_rel, 1e-8)
    opt.ftol_abs = get(prob.params, :ftol_abs, 1e-8)
    opt.xtol_rel = get(prob.params, :xtol_rel, 1e-8)
    opt.xtol_abs = get(prob.params, :xtol_abs, 1e-8)
    equality_constraint!(opt, g, get(prob.params, :tol, 0))

    if alg == :AUGLAG
        local_alg = get(prob.params, :local_alg, :LN_COBYLA)
        local_opt = Opt(local_alg, prob.params.length)
        local_opt.ftol_rel = get(prob.params, :ftol_rel, 1e-3)
        local_opt.ftol_abs = get(prob.params, :ftol_abs, 1e-3)
        local_opt.xtol_rel = get(prob.params, :xtol_rel, 1e-3)
        local_opt.xtol_abs = get(prob.params, :xtol_abs, 1e-3)
        opt.local_optimizer = local_opt
    end
    return opt
end

function SciMLBase.solve(prob::NLOptProb; MaxTime, initials)
    opt = _default_nl_opt(prob)
    opt.maxtime = MaxTime
    (optf, sol, ret) = NLopt.optimize(opt, initials)
    optsol = all_info(prob.eigfunc(sol))
    gap_der = get_gap_derivatives(prob.eigfunc, sol)
    return (; sol, optsol, gap_der...)
end

function __solve(f, alg; iterations, MaxTime=5, maxiters=1000, initials, ranges, verbosity=1, kwargs...)
    maxtime = MaxTime / iterations
    lb = map(first, ranges)
    ub = map(last, ranges)
    newinitials = map(clamp, initials, lb, ub)
    if verbosity > 0
        println("Finding sweet spot with ", alg)
        println("Initial point: ", newinitials)
    end
    prob = OptimizationProblem(f, newinitials, (1,); lb, ub)
    sol = solve(prob, alg; maxiters, maxtime, kwargs...)
    for n in 2:iterations
        newinitials = map(clamp, sol.u, lb, ub)
        verbosity > 0 && println("$n, Sweet spot:", newinitials)
        prob = OptimizationProblem(f, newinitials, (n,); lb=map(first, ranges), ub=map(last, ranges))
        sol = solve(prob, alg; maxiters, maxtime, kwargs...)
    end
    return collect(sol)
end


cost_function(gap, excgap, reduced::Number; exp, minexcgap) = cost_reduced(reduced) + cost_energy(gap, excgap; exp, minexcgap)
cost_energy(gap, excgap; minexcgap, exp) = cost_gap(gap, exp) + cost_excgap(excgap, minexcgap, exp)
cost_excgap(excgap, minexcgap, exp) = ((excgap - minexcgap) < 0 ? 1.0 + 10.0^exp * abs(excgap - minexcgap) : 0.0)

cost_gap(gap, exp) = 10.0^(exp) * abs(gap)
cost_reduced(reduced) = reduced^2

best_algs() = [BBO_probabilistic_descent(), BBO_generating_set_search(), BBO_adaptive_de_rand_1_bin_radiuslimited()]
best_alg_names() = string.([:BBO_probabilistic_descent, :BBO_generating_set_search, :BBO_adaptive_de_rand_1_bin_radiuslimited])


get_cache(c::FermionBasis, ham) = blockdiagonal(ham, c)
get_cache(::FermionBdGBasis, ham) = (Matrix(ham))

# function opt_func_cons(hamfunc, basis, ad, target; param_cost=(x...) -> 0)
#     fullsolve2(x) = fullsolve(hamfunc(x), basis)
#     cost(x, p) = target(fullsolve2(x)) + param_cost(x)
#     function cons(res, x, minexcgap)
#         sol = fullsolve2(x)
#         res[1] = sol.gap
#         res[2] = min(sol.excgap - minexcgap, 0.0)
#     end
#     if isnothing(ad)
#         return OptimizationFunction(cost; cons), fullsolve2
#     else
#         return OptimizationFunction(cost, ad; cons), fullsolve2
#     end
# end


# function SciMLBase.solve(prob::OptProb, alg::BestOf; minexcgap, initials=get_initials(prob), ranges=get_ranges(prob), exps, kwargs...)
#     fs(x) = fullsolve(prob.hamfunc(x), prob.basis)
#     function f(x, (exp, minexcgap))
#         sol = fs(x)
#         cost_function(sol.gap, sol.excgap, prob.target(sol); exp, minexcgap)
#     end
#     of = OptimizationFunction(f, AutoFiniteDiff())
#     _sols = map(alg -> _solve(of, alg; minexcgap, ranges, initials, exps, kwargs...), alg.optimizers)
#     _res = map(sol -> (; sol, optsol=fs(sol)), _sols)
#     res = filter(x -> x.optsol.excgap >= minexcgap && abs(x.optsol.gap) < 2 * 10.0^(-last(exps)), _res)
#     res = sort(res, by=x -> prob.target(x.optsol))
#     if length(res) == 0
#         @warn "No valid solutions found"
#         res = sort(_res, by=x -> prob.target(x.optsol))
#     end
#     return merge(res[1], (; all_ss=res))
# end
# function SciMLBase.solve(prob::OptProb, alg; initials=get_initials(prob), ranges=get_ranges(prob), kwargs...)
#     fs(x) = fullsolve(prob.hamfunc(x), prob.basis)
#     function f(x, (exp, minexcgap))
#         sol = fs(x)
#         cost_function(sol.gap, sol.excgap, prob.target(sol); exp, minexcgap)
#     end
#     of = OptimizationFunction(f, AutoFiniteDiff())
#     sol = _solve(of, alg; initials, ranges, kwargs...)
#     optsol = fs(sol)
#     return (; sol, optsol)
# end
# function _solve(f, alg; MaxTime=5, minexcgap=1 / 4, exps, maxiters=1000, initials, ranges, kwargs...)
#     refinements = length(exps)
#     maxtime = MaxTime / refinements
#     lb = map(first, ranges)
#     ub = map(last, ranges)
#     newinitials = map(clamp, initials, lb, ub)
#     println("Finding sweet spot with ", alg)
#     println("Initial point: ", newinitials)
#     prob = OptimizationProblem(f, newinitials, (first(exps), minexcgap); lb, ub)
#     sol = solve(prob, alg; maxiters, maxtime, kwargs...)
#     for (n, exp) in enumerate(Iterators.drop(exps, 1))
#         newinitials = map(clamp, sol.u, lb, ub)
#         println("$n, Sweet spot:", newinitials)
#         prob = OptimizationProblem(f, newinitials, (exp, minexcgap); lb=map(first, ranges), ub=map(last, ranges))
#         sol = solve(prob, alg; maxiters, maxtime, kwargs...)
#     end
#     return sol
# end

# function SciMLBase.solve(prob::OptProb, alg::IPNewton; MaxTime=5, minexcgap=1 / 4, maxiters=1000, initials=get_initials(prob), exps, kwargs...)
#     ranges = get_ranges(prob)
#     lb = map(first, ranges)
#     ub = map(last, ranges)
#     maxtime = MaxTime / (1 + final_NM)
#     println("Finding sweet spot with ", alg)
#     println("Initial point: ", initials)
#     lcons = [0.0, 0.0]
#     ucons = [0.0, Inf]
#     f, fs = opt_func(prob, alg)
#     prob2 = OptimizationProblem(f, initials, minexcgap; lb, ub, lcons, ucons, kwargs...)#, allow_f_increases = true, store_trace = true)
#     sol = solve(prob2, alg; maxiters, maxtime, kwargs...)
#     optsol = fs(sol)
#     # params = NamedTuple(zip((:Δ, :δϕ, :ε), decompose(sol)))
#     return (; alg, sol, optsol)
#     # return sol
# end