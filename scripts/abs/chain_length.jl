using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra, BlackBoxOptim
using Plots
using Symbolics

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
fixedparams = (; t=0.5, θ=parameter(π / 2, :diff), ϕ=0, V=0, Δ=0, U=1 * 10, h=2)
c = FermionBdGBasis(1:N, (:↑, :↓))
c = FermionBasis(1:N, (:↑, :↓), qn=QuantumDots.parity)
@variables Δ1::Real, μs[1:div(N + 1, 2)]::Real
@time symbolic_ham = hamiltonian(c; Δ1, μ=reflect(μs, N), fixedparams...);
@which build_function(symbolic_ham, [Δ1, μs...], expression=Val{false})
# cache = hamfuncs[1]([1, ones(div(N+1,2)...])
##
function get_hamfuncs(c, fixedparams)
    N = div(QuantumDots.nbr_of_fermions(c), 2)
    @variables Δ1::Real, μs[1:div(N + 1, 2)]::Real
    symbolic_ham = hamiltonian(c; Δ1, μ=reflect(μs, N), fixedparams...)
    hamfuncs = build_function(symbolic_ham, [Δ1, μs...], expression=Val{false})
    ps = rand(1 + div(N + 1, 2))
    cache = hamfuncs[1](ps)
    hamfuncs[2](cache, ps)
    return hamfuncs, cache
end
@code_warntype get_hamfuncs(FermionBasis(1:4, (:↑, :↓), qn=QuantumDots.parity), fixedparams);

# (h1,h2),cache = get_hamfuncs(FermionBdGBasis(1:2, (:↑, :↓)), fixedparams)
# h1(cache, [1,2]) |> QuantumDots.BdGMatrix
# h1(cache, [1,2])

##
fixedparams = (; t=0.5, θ=parameter(π / 2.0, :diff), ϕ=0.0, V=0.0, Δ=0.0, U=10.0, h=2.0)
function optimize(N, fixedparams; basis=:bdg, MaxTime=5, minexcgap=1 / 4, exps=collect(range(0.1, 3, length=4)), TargetFitness=1e-6, initials=[1.0, (fixedparams.h .* ones(div(N + 1, 2)))...], tracemode=:silent, transport=missing)
    c = basis == :bdg ? FermionBdGBasis(1:N, (:↑, :↓)) : FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
    hmax = 3.0 * (fixedparams.h + fixedparams.U)
    hamfuncs, cache = get_hamfuncs(c, fixedparams)
    target = x -> MPU(x) + LDf(x)
    opt = LongerPoorMansMajoranas.Optimizer(;
        hamfunc=(Δ1, μs...) -> (hamfuncs[2](cache, [Δ1, reflect(μs, N)...]); return cache),
        ranges=[(-1.0, 1.0), [(-hmax / 2, hmax) for i in 1:div(N + 1, 2)]...],
        initials,
        MaxTime = Float64(MaxTime), minexcgap,
        exps, TargetFitness,
        tracemode,
        target,
        basis=c)

    ss = get_sweet_spot(opt)
    bc = best_candidate(ss)
    optsol = fullsolve(opt.hamfunc(bc...), c; transport)
    params = NamedTuple(zip((:Δ1, :μ), (first(bc), bc[2:end])))
    (; ss, bc, optsol, params)
end
##
optsols = [optimize(N, fixedparams; MaxTime=2, tracemode = :compact) for N in 2:4];
optsols2 = [optimize(N, fixedparams; MaxTime=N^2, basis=:mb, initials=sol.bc,tracemode = :compact) for (N, sol) in zip(2:4, optsols)];
optsols3 = [optimize(N, fixedparams; MaxTime=N^2, basis=:mb, initials=sol.bc) for (N, sol) in zip(2:4, optsols)]
##
sols = optsols2
for sol in sols
    plotmajcoeffs(sol.optsol.majcoeffs...) |> display
end
map(sol -> LDf(sol.optsol), sols)
map(sol -> MPU(sol.optsol), sols)
map(sol -> sol.bc, sols)
map(sol -> sol.optsol.excgap, sols)
##
for sol in optsols2
    plotmajcoeffs(sol.optsol.majcoeffs...) |> display
end
map(sol -> LDf(sol.optsol), optsols2)
map(sol -> MPU(sol.optsol), optsols2)
map(sol -> sol.bc, optsols2)
map(sol -> sol.optsol.excgap, optsols2)
