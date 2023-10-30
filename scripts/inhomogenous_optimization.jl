using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra, BlackBoxOptim
using Plots
##
N = 2
c1 = FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
c2 = FermionBdGBasis(1:N, (:↑, :↓))
transport = Transport(QuantumDots.Pauli(), (; T=1 / 40, μ=(0.0, 0.0)))
##
fixedparams = (; t=0.5, θ=parameter(π / 2, :diff), ϕ=0, V=0, Δ=1, U=0 * 10, h=2)
c = c2
opt = LongerPoorMansMajoranas.Optimizer(
    hamfunc=(Δ1, μs...) -> hamiltonian(c; Δ1, μ=μs, fixedparams...),
    ranges=[(0.01, 1.0), [(-1 * fixedparams.h, 3 * fixedparams.h) for i in 1:N]...],
    initials=[1.0, (fixedparams.h .+ zeros(N))...];
    MaxTime=2, minexcgap=1 / 4,
    exps=collect(range(0.1, 3, length=4)),
    tracemode=:compact,
    target=LD,
    basis=c)

ss = get_sweet_spot(opt)
bc = best_candidate(ss)
optsol = fullsolve(opt.hamfunc(bc...), c; transport)

csdata1 = charge_stability_scan(merge(fixedparams, NamedTuple(zip((:Δ1, :μ), (first(bc), bc[2:end])))), 0.1*5fixedparams.h, 0.1*5fixedparams.h; basis=c1)
csdata2 = charge_stability_scan(merge(fixedparams, NamedTuple(zip((:Δ1, :μ), (first(bc), bc[2:end])))), 0.1*5fixedparams.h, 0.1*5fixedparams.h; basis=c2)
##
drho = optsol.reduced.cell_matrices.even[1] -
optsol.reduced.cell_matrices.odd[1]

##
csdata = csdata2
heatmap(csdata[:μ1], csdata[:μ2], map(x -> tanh(x.gap), csdata[:data]); clims=(-1, 1), c=:balance)
heatmap(csdata[:μ1], csdata[:μ2], map(MPU, csdata[:data]); clims=(0, 1), c=:viridis)
heatmap(csdata[:μ1], csdata[:μ2], map(LDf, csdata[:data]); clims=(0, 1), c=:viridis)
heatmap(csdata[:μ1], csdata[:μ2], map(LD, csdata[:data]); clims=(0, 2), c=:viridis)


##
[map(LD, csdata1[:data]) - map(LD, csdata2[:data])x |> norm for x in range(.75,.76,1000)] |> findmin
map(MP, csdata1[:data]) - map(MP, csdata2[:data]) |> norm
map(LDf, csdata1[:data]) - map(LDf, csdata2[:data]) |> norm