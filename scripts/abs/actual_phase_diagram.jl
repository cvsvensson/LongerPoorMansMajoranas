using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra#, BlackBoxOptim
using Plots
using Symbolics
using SkewLinearAlgebra

##
function hoppingH(k, ε, (t1, t2))
    t = 2 * (cospi(k) * t1 + cospi(2k) * t2)
    [ε+t 0; 0 ε+t]
end
function pairingH(k, (Δ1, Δ2))
    Δ = 2 * 1im * (sinpi(k) * Δ1 + sinpi(2k) * Δ2)
    [0 Δ; -Δ 0]
end
bdgH(k, ε, (t1, t2), (Δ1, Δ2); kwargs...) = BdGMatrix(hoppingH(k, ε, (t1, t2)), pairingH(k, (Δ1, Δ2)); kwargs...)[2:3, 2:3]
skewH(k, ε, (t1, t2), (Δ1, Δ2); kwargs...) = QuantumDots.bdg_to_skew(bdgH(k, ε, (t1, t2), (Δ1, Δ2); kwargs...))
function topoQ(ε, (t1, t2), (Δ1, Δ2); kwargs...)
    pf1 = pfaffian(skewH(0, ε, (t1, t2), (Δ1, Δ2); kwargs...))
    pf2 = pfaffian(skewH(1, ε, (t1, t2), (Δ1, Δ2); kwargs...))
    sign(pf1 * pf2)
end
##

h = bdgH(1, 0, (1, 0), (1, 2))
skewH(1, 0, (1, 0), (1, 2))
pfaffian(skewH(0, 0, (1, 0), (1, 0)))

##

ϵs = range(-4, 4, length=100)
ts = range(0, 2, length=101)
ks = range(-1, 1, length=100)
data = [topoQ(eps, (t, .2), (1, 0)) for eps in ϵs, t in ts]
heatmap(ts, ϵs, data, clims=(-1, 1), c=:redsblues, xlabel="t", ylabel="ε", title="Topological invariant")
##
let t = 1
    energies = [[eigvals(bdgH(k, eps, (1, 0.0), (1.1, 0.0))) for k in ks] for eps in [-2t, 0, 2t]]
    plot([plot(ks, stack(es)') for es in energies]...)
end
##
bdgH(0.3, -1, (0.5, 0.0), (1, 0))
skewH(0, 0.2, (0.5, 0.0), (1, 0))
