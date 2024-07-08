using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra#, BlackBoxOptim
using Plots
using Symbolics
using SkewLinearAlgebra
using StaticArrays
using GellMannMatrices
using Roots
using ForwardDiff
##
function getHs(k, ε, (t1, t2), (Δ1, Δ2))
    s1, c1 = sincos(k)
    s2, c2 = sincos(2k)
    h0 = -imag(t1) * s1 - imag(t2) * s2
    h1 = -imag(Δ1) * s1 - imag(Δ2) * s2
    h2 = -real(Δ1) * s1 - real(Δ2) * s2
    h3 = real(t1) * c1 + real(t2) * c2 + ε / 2
    @SVector [h0, h1, h2, h3]
end
function bdgH(k, ε, (t1, t2), (Δ1, Δ2))#, paulis=paulis)
    hs = getHs(k, ε, (t1, t2), (Δ1, Δ2))
    Hermitian(mapreduce(*, +, hs, paulis))
end
function bdgQ(k, ε, (t1, t2), (Δ1, Δ2), paulis=paulis)
    hs = getHs(k, ε, (t1, t2), (Δ1, Δ2))
    @SMatrix [hs[1]+hs[2] hs[4]+1im*hs[3]; hs[4]-1im*hs[3] hs[1]-hs[2]]
end
@inline function fast_energies(m::Hermitian{<:Any,SMatrix{2,2,T,L}}) where {T,L}
    t = tr(m)
    d = m[1] * m[4] - m[2] * m[3] #det(m)
    s = sqrt(t^2 - 4d)
    (real(t + s) / 2, real(t - s) / 2)
end
##
@variables k t::Real Δ::Real Ez::Real θ::Real t1::Complex t2::Complex Δ1::Complex Δ2::Complex ε::Real δϕ::Real δε::Real
##
const paulis = SVector{4}(pushfirst!(map(SMatrix{2,2}, gellmann(2)), I(2)))
h = bdgH(k, ε, (t1, t2), (Δ1, Δ2))
[tr(h * σ) for σ in paulis]
h0 = substitute.(h, k => 0)
hpi = substitute.(h, k => pi)
[tr(h0 * σ) for σ in paulis]
[tr(hpi * σ) for σ in paulis]
##
skewH(k, ε, (t1, t2), (Δ1, Δ2); kwargs...) = QuantumDots.bdg_to_skew(bdgH(k, ε, (t1, t2), (Δ1, Δ2)); kwargs...)
substitute.(skewH(k, ε, (t1, t2), (Δ1, Δ2); check=false), k => 0)
function topoQ(ε, (t1, t2), (Δ1, Δ2); kwargs...)
    sign((2real(t1 + t2) + ε) * (2real(-t1 + t2) + ε))
end
function energy_gap(ε, (t1, t2), (Δ1, Δ2))
    f(k) = fast_energies(bdgH(k, ε, (t1, t2), (Δ1, Δ2)))[1]
    Df(k) = ForwardDiff.derivative(f, k)
    E0 = min(abs(f(0)), abs(f(pi)))
    sol1 = find_zeros(f, -pi, pi)
    if length(sol1) > 0
        return min(minimum(abs ∘ f, sol1), E0)
    end
    sol2 = find_zeros(Df, -pi, pi)
    return min(minimum(abs ∘ f, sol2), E0)
end
##
@time energy_gap(0.5, (2exp(1im * pi / 4), 0), (1, 0))
@code_warntype eigvals(bdgH(0.5, 0.5, (2exp(1im * pi / 4), 0), (1, 0)))[1]
##
let ks = range(-pi, pi + 0.1, 100), f(k) = collect(eigvals(bdgH(k, 2, (exp(0.0 * 1im * pi / 4), 0), (1, 0)))), Df
    Df = x -> ForwardDiff.derivative(f, x)
    plot(plot(ks, stack([f(k) for k in ks])'),
        plot(ks, stack([Df(k) for k in ks])'))
end
##
@time [energy_gap(eps, (exp(1im * 1) * t, 0), (1, 0)) for eps in range(-4, 4, length=10), t in range(0, 4, length=11)];
@profview [energy_gap(eps, (exp(1im * 1) * t, 0), (1, 0)) for eps in range(-4, 4, length=100), t in range(0, 4, length=101)];

## Standard kitaev phase diagram
let ϕ = 0.0 * pi / 10
    ϵs = range(-4, 4, length=100)
    ts = range(0, 4, length=101)
    Δ = (1, 0)
    @time data = [topoQ(eps, (exp(1im * ϕ) * t, 0), Δ; check=false) for eps in ϵs, t in ts]
    @time data2 = [energy_gap(eps, (exp(1im * ϕ) * t, 0), Δ) for eps in ϵs, t in ts]
    prob = BPProblem((k, p) -> Matrix(bdgH(k, p[2], (exp(1im * ϕ) * p[1], 0), Δ)))
    @time calcPhaseDiagram(prob, ts, ϵs; plot=true)
    p1 = heatmap(ts, ϵs, data, clims=(-1, 1), c=:redsblues, xlabel="t", ylabel="ε", title="Topological invariant")
    p2 = heatmap(ts, ϵs, data2, clims=(0, 1), c=:viridis, xlabel="t", ylabel="ε", title="Energy gap")
    plot(p1, p2)
end
## Match 1303.3304
let t = 1, ϕ = 0.5 * pi / 2
    ϵs = range(-4, 4, length=100)
    Δs = range(-2, 2, length=101)
    data = [topoQ(eps, (exp(1im * ϕ) * t, 0), (Δ, 0); check=false) for eps in ϵs, Δ in Δs]
    data2 = [energy_gap(eps, (exp(1im * ϕ) * t, 0), (Δ, 0)) for eps in ϵs, Δ in Δs]
    p1 = heatmap(ϵs, Δs, data', clims=(-1, 1), c=:redsblues, xlabel="ε", ylabel="Δ", title="Topological invariant")
    p2 = heatmap(ϵs, Δs, data2', clims=(0, 1), c=:viridis, xlabel="ε", ylabel="Δ", title="Energy gap")
    plot(p1, p2)
end
##
let t = 1
    ks = range(-pi, pi, length=100)
    energies = [[fast_energies(bdgH(k, eps, (1, 0.0), (1.1, 0.0))) for k in ks] for eps in [-2t, 0, 2t]]
    plot([plot(ks, stack(es)') for es in energies]...)
end
##
let t = 1
    ks = range(-pi, pi, length=100)
    energies = [[fast_energies(bdgH(k, eps, (1im, 0.0), (1.1, 0.0))) for k in ks] for eps in [-2t, 0, 2t]]
    plot([plot(ks, stack(es)') for es in energies]...)
end
##
function effective_coeffs(; ε, t, θ, Δ, Ez, δϕ)
    Δ_aa, Δ_bb, Δ_ab, Δ_ba, t_aa, t_ab, t_ba, t_bb = LongerPoorMansMajoranas.perturbative_coeffs_homogeneous(; ε, t, θ, Δ, δϕ)
    ε0 = Ez - sqrt(Δ^2 + ε^2)
    (; E, ε_ba, ε_ab, t_nn, Δ_nn) = LongerPoorMansMajoranas.second_order_coeffs_homogeneous(; Δ_ba, Δ_ab, t_ba, t_ab, Ez, ε, Δ)
    (; ε=real(ε_ba + ε_ab) + ε0, t1=t_aa, t2=t_nn, Δ1=-conj(Δ_aa), Δ2=Δ_nn)
end
function topoQ(; first_order=false, ε, t, θ, Δ, Ez, δϕ)
    (; ε, t1, t2, Δ1, Δ2) = effective_coeffs(; ε, t, θ, Δ, Ez, δϕ)
    topoQ(ε, first_order ? (t1, 0t2) : (t1, t2), first_order ? (Δ1, 0Δ2) : (Δ1, Δ2); check=false)
end
function energy_gap(; first_order=false, ε, t, θ, Δ, Ez, δϕ, kwargs...)
    (; ε, t1, t2, Δ1, Δ2) = effective_coeffs(; ε, t, θ, Δ, Ez, δϕ)
    energy_gap(ε, first_order ? (t1, 0t2) : (t1, t2), first_order ? (Δ1, 0Δ2) : (Δ1, Δ2); kwargs...)
end
function bdgH(k; first_order=false, ε, t, θ, Δ, Ez, δϕ)
    (; ε, t1, t2, Δ1, Δ2) = effective_coeffs(; ε, t, θ, Δ, Ez, δϕ)
    bdgH(k, ε, first_order ? (t1, 0t2) : (t1, t2), first_order ? (Δ1, 0Δ2) : (Δ1, Δ2))
end
##
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), Δ=1, Ez=3)
ϵs = fixedparams.Ez - 0.1 .+ 0.4 .* range(-1, 1, length=150)
δϕs = range(0, pi, length=151)
##
dataE1 = [energy_gap(; ε, first_order=true, δϕ, fixedparams...) for ε in ϵs, δϕ in δϕs]
dataQ1 = [topoQ(; first_order=true, ε, δϕ, fixedparams...) for ε in ϵs, δϕ in δϕs]
dataE2 = [energy_gap(; ε, first_order=false, δϕ, fixedparams...) for ε in ϵs, δϕ in δϕs]
dataQ2 = [topoQ(; first_order=false, ε, δϕ, fixedparams...) for ε in ϵs, δϕ in δϕs]
prob1 = BPProblem((k, p) -> Matrix(bdgH(k; first_order=true, ε=p[2], δϕ=p[1], fixedparams...)))
prob2 = BPProblem((k, p) -> Matrix(bdgH(k; first_order=false, ε=p[2], δϕ=p[1], fixedparams...)))
@time dataTQ1 = calcPhaseDiagram(prob1, δϕs, ϵs; plot=true)
@time dataTQ2 = calcPhaseDiagram(prob2, δϕs, ϵs; plot=true)
##
p1E = heatmap(ϵs, δϕs, dataE1', clims=(1e-10, 1e-2), colorbar_scale=:identity, c=:viridis, xlabel="ε", ylabel="δϕ", title="Energy gap 1st order")
p1Q = heatmap(ϵs, δϕs, dataQ1', clims=(-1, 1), c=:redsblues, xlabel="ε", ylabel="δϕ", title="Topological invariant")
p2E = heatmap(ϵs, δϕs, dataE2', clims=(1e-10, 1e-2), colorbar_scale=:identity, c=:viridis, xlabel="ε", ylabel="δϕ", title="Energy gap 2nd order")
p2Q = heatmap(ϵs, δϕs, dataQ2', clims=(-1, 1), c=:redsblues, xlabel="ε", ylabel="δϕ", title="Topological invariant")
plot(p1E, p1Q, p2E, p2Q)
##
a = FermionBdGBasis(1:5)
c = FermionBdGBasis(1:5, (:↑, :↓))
f, f!, cache = hamfunc(Hδϕ_Hε(), c, merge(fixedparams, (; U=0, V=0)))
data_p1 = [(
    begin
        H = LongerPoorMansMajoranas.perturbative_hamiltonian_homogeneous(a, 1; ε, δϕ, fixedparams...)
        fullsolve(H, a)
    end
) for ε in ϵs, δϕ in δϕs]
data_p2 = [(
    begin
        H = LongerPoorMansMajoranas.perturbative_hamiltonian_homogeneous(a, 2; ε, δϕ, fixedparams...)
        fullsolve(H, a)
    end
) for ε in ϵs, δϕ in δϕs]
data_f = [(
    begin
        H = f!(cache, [δϕ, ε])
        fullsolve(H, c)
    end
) for ε in ϵs, δϕ in δϕs]
##
p3 = heatmap(ϵs, δϕs, map(LDbdg, data_p2)', c=:viridis, xlabel="ε", ylabel="δϕ", title="LD, 2nd order, 20 sites", clims=(0, 1))
p4 = heatmap(ϵs, δϕs, map(LDbdg, data_f)', c=:viridis, xlabel="ε", ylabel="δϕ", title="LD, full model, 20 sites", clims=(0, 1))
p2nd = plot(p2E, p2Q, p3, p4, layout=(2, 2), size=(800, 500))
##
p3 = heatmap(ϵs, δϕs, map(LDbdg, data_p1)', c=:viridis, xlabel="ε", ylabel="δϕ", title="LD, 1st order, 20 sites")
p4 = heatmap(ϵs, δϕs, map(LDbdg, data_f)', c=:viridis, xlabel="ε", ylabel="δϕ", title="LD, full model, 20 sites")
p1st = plot(p1E, p1Q, p3, p4, layout=(2, 2), size=(800, 500))
##
savefig(p2nd, "phase_diagrams_perturbative_vs_full_20.png")
##
heatmap([(
        begin
            (; ε, t1, t2, Δ1, Δ2) = effective_coeffs(; fixedparams..., ε, δϕ)
            real(t1 - Δ1)
        end
    ) for ε in ϵs, δϕ in δϕs]', c=:redsblues, clims=0.01 .* (-1, 1))
##

effective_coeffs(; ε, fixedparams..., δϕ)
let e = ε
    (; ε, t1, t2, Δ1, Δ2) = map(simplify, effective_coeffs(; ε=e, Ez, θ=parameter(θ, :diff), t, Δ, δϕ))
    substitute(map(simplify, bdgQ(k, ε, (t1, t2), (Δ1, Δ2)))[1, 2], k => 0)
end

##
f = plot()
for s in [:t1, :t2, :Δ1, :Δ2]
    x, y = eachrow(map(x -> x.val, stack(reim.([substitute(coeffs[s], δϕ => x) for x in range(0, sqrt(pi), 20) .^ 2]))))
    scatter!(f, x, y, aspectratio=1, label=string(s))
end
f
##
f = plot()
for s in [:t1, :Δ1, :t2, :Δ2]
    ϕs = range(0, pi, 20)
    D = Differential(δϕ)
    # expression = expand_derivatives(D(coeffs[s]))
    expression = expand_derivatives(angle(D(coeffs[s]))) |> simplify
    println(expression)
    # a = map(x -> x.val, (angle.([substitute(expression, δϕ => x) for x in ϕs])))
    a = map(x -> x.val + 0.1rand(), (([substitute(expression, δϕ => x) for x in ϕs])))
    plot!(f, ϕs, a, markers=true, label=string(s))
end
f

##

expression = effective_coeffs(; ε, Ez, θ=parameter(θ, :diff), t, Δ, δϕ) |> simplify
D = Differential(δϕ)
for s in [:t1, :Δ1]#, :t2, :Δ2]
    expression2 = expand_derivatives(D(expression[s])) |> simplify
    expression3 = expand_derivatives(D(angle(expression2))) |> simplify
    println("Derivative for $s:", expression3)
end
##
for s in [:t1, :Δ1]#, :t2, :Δ2]
    expression2 = expand_derivatives(D(expression[s])) |> simplify
    Symbolics.coeff(expression2, cos(δϕ))
    # expression3 = expand_derivatives(D(angle(expression2))) |> simplify
    println("Derivative for $s:", expression3)
end

##
expression2 = expand_derivatives(D(expression[:t1])) |> simplify
nn, dd = Symbolics.arguments(Symbolics.value(real(expression2)))
Symbolics.coeff(simplify(nn), cos(δϕ))