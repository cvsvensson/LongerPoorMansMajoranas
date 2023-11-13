using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra, BlackBoxOptim
using Plots
function plotmajcoeffs(ws, zs)
    k = string.(keys(ws).values)
    labels = ["w", "z"]
    plot(k, [abs2.(ws).values, abs2.(zs).values], label=labels, legend=:topleft)
end

##
N = 2
c1 = FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
c2 = FermionBdGBasis(1:N, (:↑, :↓))
transport = Transport(QuantumDots.Pauli(), (; T=1 / 40, μ=(0.0, 0.0)))
target = x -> MPU(x) + LDf(x)
##
fixedparams = (; t=0.5, θ=parameter(π / 2, :diff), ϕ=0, V=0, Δ=0, U=0 * 10, h=2)
c = c2
opt = LongerPoorMansMajoranas.Optimizer(;
    hamfunc=(Δ1, μs...) -> hamiltonian(c; Δ1, μ=reflect(μs, N), fixedparams...),
    ranges=[(-1.0, 1.0), [(-1 * fixedparams.h, 3 * fixedparams.h) for i in 1:div(N + 1, 2)]...],
    initials=[1.0, (fixedparams.h .* ones(div(N + 1, 2)))...],
    MaxTime=5, minexcgap=1 / 4,
    exps=collect(range(0.1, 3, length=4)),
    tracemode=:compact,
    target,
    basis=c)

ss = get_sweet_spot(opt)
bc = best_candidate(ss)
optsol = fullsolve(opt.hamfunc(bc...), c; transport)
plotmajcoeffs(optsol.majcoeffs...)

##
csdata1 = charge_stability_scan(merge(fixedparams, NamedTuple(zip((:Δ1, :μ), (first(bc), bc[2:end])))), 5fixedparams.h, 5fixedparams.h; basis=c1)
csdata2 = charge_stability_scan(merge(fixedparams, NamedTuple(zip((:Δ1, :μ), (first(bc), bc[2:end])))), 5fixedparams.h, 5fixedparams.h; basis=c2)
##
drho = optsol.reduced.cell_matrices.even[1] -
       optsol.reduced.cell_matrices.odd[1]

##
csdata = csdata2
heatmap(csdata[:μ1], csdata[:μ2], map(x -> tanh(x.gap), csdata[:data]); clims=(-1, 1), c=:balance)
heatmap(csdata[:μ1], csdata[:μ2], map(MPU, csdata[:data]); clims=(0, 1), c=:viridis)
heatmap(csdata[:μ1], csdata[:μ2], map(tanh ∘ LDf, csdata[:data]); clims=(0, 1), c=:viridis)
heatmap(csdata[:μ1], csdata[:μ2], map(tanh ∘ LD, csdata[:data]); clims=(0, 1), c=:viridis)
heatmap(csdata[:μ1], csdata[:μ2], map(tanh ∘ target, csdata[:data]); clims=(0, 1), c=:viridis)


##
map(MP, csdata1[:data]) - map(MP, csdata2[:data]) |> norm #Should be close to zero without interactions
map(LDf, csdata1[:data]) - map(LDf, csdata2[:data]) |> norm #Should be close to zero without interactions
map(LD, csdata1[:data]) - map(LD, csdata2[:data]) |> norm #May not be zero

## Inhomogeneous Δ1 
# fixedparams = (; t=0.5, θ=parameter(π / 2, :diff), ϕ=0, V=0, Δ=0, U=0 * 10, h=2)
c = c2
Nhalf = div(N + 1, 2)
# initials=[ones(N - 1)..., (fixedparams.h .+ zeros(Nhalf))...]
initials = [[bc[1] for i in 1:div(N, 2)]..., bc[2:end]...]
opt = LongerPoorMansMajoranas.Optimizer(;
    hamfunc=(Δ1μs...) -> hamiltonian(c; Δ1=reflect(Δ1μs[1:N÷2], N - 1; pad=[0.0]), μ=reflect(Δ1μs[(1+N÷2):end], N), fixedparams...),
    ranges=[[(-1.0, 1.0) for i in 1:N÷2]..., [(-1 * fixedparams.h, 3 * fixedparams.h) for i in 1:Nhalf]...],
    initials,
    MaxTime=20, minexcgap=1 / 4,
    exps=collect(range(.1, 3, length=3)),
    tracemode=:compact,
    # target=x -> MPU(x) + LDf(x)^2,
    target,
    basis=c)

ss2 = get_sweet_spot(opt)
bc2 = best_candidate(ss2)
optsol2 = fullsolve(opt.hamfunc(bc2...), c; transport)
plotmajcoeffs(optsol2.majcoeffs...)
##
newparams = merge(fixedparams, NamedTuple(zip((:Δ1, :μ), (reflect(bc2[1:N÷2], N - 1; pad=[0.0]), reflect(bc2[N÷2+1:end], N)))))
csdata12 = charge_stability_scan(newparams, 5fixedparams.h, 5fixedparams.h, 50; basis=c1, transport)
csdata22 = charge_stability_scan(newparams, 5fixedparams.h, 5fixedparams.h; basis=c2)
##
map(MP, csdata12[:data]) - map(MP, csdata22[:data]) |> norm #Should be close to zero without interactions
map(LDf, csdata12[:data]) - map(LDf, csdata22[:data]) |> norm #Should be close to zero without interactions
map(LD, csdata12[:data]) - map(LD, csdata22[:data]) |> norm #May not be zero

##
csdata = csdata12
heatmap(csdata[:μ1], csdata[:μ2], map(x -> tanh(x.gap), csdata[:data]); clims=(-1, 1), c=:balance)
heatmap(csdata[:μ1], csdata[:μ2], map(MPU, csdata[:data]); clims=(0, 1), c=:viridis)
heatmap(csdata[:μ1], csdata[:μ2], map(MP, csdata[:data]); clims=(0, 1), c=:viridis)
heatmap(csdata[:μ1], csdata[:μ2], map(tanh ∘ LDf, csdata[:data]); clims=(0, 1), c=:viridis)
heatmap(csdata[:μ1], csdata[:μ2], map(tanh ∘ LD, csdata[:data]); clims=(0, 1), c=:viridis)
heatmap(csdata[:μ1], csdata[:μ2], map(tanh ∘ target, csdata[:data]); clims=(0, 1), c=:viridis)

heatmap(csdata[:μ1], csdata[:μ2], map(x -> x.parity, csdata[:data]); clims=(-1, 1), c=:berlin)
##
heatmap(csdata[:μ1], csdata[:μ2], map(x -> tanh(real(x.conductance[1, 1])), csdata[:data]); clims=(0, 1), c=:amp)
heatmap(csdata[:μ1], csdata[:μ2], map(x -> tanh(real(x.conductance[1, 2])), csdata[:data]); clims=(-1, 1), c=:balance)

