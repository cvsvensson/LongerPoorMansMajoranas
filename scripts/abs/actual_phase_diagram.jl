using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra#, BlackBoxOptim
using Plots
using Symbolics
using SkewLinearAlgebra
using Optimization, OptimizationBBO
##
function hoppingH(k, ε, (t1, t2))
    t = 2 * real(cis(k) * t1 + cis(2k) * t2)
    [ε+t 0; 0 ε+t]
end
function pairingH(k, (Δ1, Δ2))
    Δ = 2 * 1im * imag(cis(k) * Δ1 + cis(2k) * Δ2)
    [0 Δ; -Δ 0]
end
bdgH(k, ε, (t1, t2), (Δ1, Δ2); kwargs...) = BdGMatrix(hoppingH(k, ε, (t1, t2)), pairingH(k, (Δ1, Δ2)); kwargs...)[2:3, 2:3]
skewH(k, ε, (t1, t2), (Δ1, Δ2); kwargs...) = QuantumDots.bdg_to_skew(bdgH(k, ε, (t1, t2), (Δ1, Δ2); kwargs...); kwargs...)
function topoQ(ε, (t1, t2), (Δ1, Δ2); kwargs...)
    pf1 = pfaffian(skewH(0, ε, (t1, t2), (Δ1, Δ2); kwargs...))
    pf2 = pfaffian(skewH(1, ε, (t1, t2), (Δ1, Δ2); kwargs...))
    sign(pf1 * pf2)
end
function energy_gap(ε, t, Δ; maxtime=0.0001)
    # f(k, p) = abs(eigvals!(bdgH(only(k), ε, t, Δ), sortby=abs)[1])
    f(k, p) = sqrt(real(tr(bdgH(only(k), ε, t, Δ)^2)))
    # find_zero(f, 0.0, 2pi; maxtime)
    prob = OptimizationProblem(f, [0.0]; lb=-pi, ub=pi)
    sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxtime)
    f(sol[1], nothing)
end
##

h = bdgH(1, 0, (1 + 1im, 0), (1, 2))
skewH(1, 0, (1 + 1im, 0), (1, 2))
pfaffian(skewH(0, 0, (1, 0), (1, 0)))

##

ϵs = range(-4, 4, length=50)
ts = range(0, 2, length=51)
ks = range(-1, 1, length=100)
data = [topoQ(eps, (t, 0), (1, 0)) for eps in ϵs, t in ts]
data2 = [energy_gap(eps, (t, 0), (1, 0)) for eps in ϵs, t in ts]
heatmap(ts, ϵs, data, clims=(-1, 1), c=:redsblues, xlabel="t", ylabel="ε", title="Topological invariant")
heatmap(ts, ϵs, data2, clims=(0, 1), c=:viridis, xlabel="t", ylabel="ε", title="Energy gap")

##
let t = 1
    energies = [[eigvals(bdgH(k, eps, (1im, 0.0), (1.1, 0.0))) for k in ks] for eps in [-2t, 0, 2t]]
    plot([plot(ks, stack(es)') for es in energies]...)
end
##
let t = 1
    energies = [[eigvals(bdgH(k, eps, (1im, 0.0), (1.1, 0.0))) for k in ks] for eps in [-2t, 0, 2t]]
    plot([plot(ks, stack(es)') for es in energies]...)
end
##
bdgH(0.3, -1, (0.5, 0.0), (1, 0))
skewH(0, 0.2, (0.5, 0.0), (1, 0))

##
function topoQ(; ε, t, θ, Δ, Ez, δϕ)
    # coeffs = stack(LongerPoorMansMajoranas.perturbative_coeffs(n; Δ, ε, θ, δϕ, t) for n in 1:2)
    # Δ_aa, Δ_bb, Δ_ab, Δ_ba, t_aa, t_ab, t_ba, t_bb = collect(eachrow(coeffs))
    # (; E1, E2, ε_ba, ε_ab, t_nn, Δ_nn) = LongerPoorMansMajoranas.second_order_coeffs(1; Δ_ba, Δ_ab, t_ba, t_ab, Ez, ε, Δ)
    Δ_aa, Δ_bb, Δ_ab, Δ_ba, t_aa, t_ab, t_ba, t_bb = LongerPoorMansMajoranas.perturbative_coeffs(1; ε, t, θ, Δ, δϕ)
    (; E1, E2, ε_ba, ε_ab, t_nn, Δ_nn) = LongerPoorMansMajoranas.second_order_coeffs(1; Δ_ba=[Δ_ba, Δ_ba], Δ_ab=[Δ_ab, Δ_ab], t_ba=t_ba * [1, 1], t_ab=t_ab * [1, 1], Ez, ε, Δ)
    t1 = t_aa
    t2 = t_nn
    Δ1 = Δ_aa
    Δ2 = Δ_nn
    topoQ(real(ε_ba + ε_ab), (t1, t2), (Δ1, Δ2))
end
function energy_gap(; ε, t, θ, Δ, Ez, δϕ, maxtime=0.0001)
    Δ_aa, Δ_bb, Δ_ab, Δ_ba, t_aa, t_ab, t_ba, t_bb = LongerPoorMansMajoranas.perturbative_coeffs(1; ε, t, θ, Δ, δϕ)
    (; E1, E2, ε_ba, ε_ab, t_nn, Δ_nn) = LongerPoorMansMajoranas.second_order_coeffs(1; Δ_ba=[Δ_ba, Δ_ba], Δ_ab=[Δ_ab, Δ_ab], t_ba=t_ba * [1, 1], t_ab=t_ab * [1, 1], Ez, ε, Δ)
    t1 = t_aa
    t2 = t_nn
    Δ1 = Δ_aa
    Δ2 = Δ_nn
    energy_gap(ε_ba + ε_ab, (t1, t2), (Δ1, Δ2); maxtime)
end
topoQ(; ε=[0, 0, 0], δϕ=[0, 0, 0], fixedparams...)
##
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), Δ=[1, 1, 1], Ez=3)
ϵs = fixedparams.Ez .+ range(-5, 5, length=50)
δϕs = range(0, pi, length=51)
data = [energy_gap(; ε=eps * [1, 1, 1], δϕ=δϕ * [1, 1], fixedparams...) for eps in ϵs, δϕ in δϕs]
heatmap(ϵs, δϕs, data', clims=(0, 0.2), c=:viridis, xlabel="ε", ylabel="δϕ", title="energy gap")
##
data = [topoQ(; ε=eps * [1, 1, 1], δϕ=δϕ * [1, 1], fixedparams...) for eps in ϵs, δϕ in δϕs]
heatmap(δϕs, ϵs, data, clims=(-1, 1), c=:redsblues, xlabel="δϕ", ylabel="ε", title="Topological invariant")
##