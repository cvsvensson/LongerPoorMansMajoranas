using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra#, BlackBoxOptim
using Plots
using Symbolics
using SkewLinearAlgebra
using Optimization, OptimizationBBO
using StaticArrays
using GellMannMatrices
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
function bdgH(k, ε, (t1, t2), (Δ1, Δ2), paulis=paulis)
    hs = getHs(k, ε, (t1, t2), (Δ1, Δ2))
    mapreduce(*, +, hs, paulis)
    # h0*paulis[1] + h1*paulis[2] + h2*paulis[3] + h3*paulis[4]
    # sum(v .* paulis)
    # mapreduce(*, +, v, paulis)
end
function bdgQ(k, ε, (t1, t2), (Δ1, Δ2), paulis=paulis)
    hs = getHs(k, ε, (t1, t2), (Δ1, Δ2))
    @SMatrix [hs[1]+hs[2] hs[4]+1im*hs[3]; hs[4]-1im*hs[3] hs[1]-hs[2]]
end
##
@variables k t::Real Δ::Real Ez::Real θ::Real t1::Complex t2::Complex Δ1::Complex Δ2::Complex ε::Real δϕ::Real
##
const paulis = SVector{4}(pushfirst!(map(SMatrix{2,2}, gellmann(2)), I(2)))
h = bdgH(k, ε, (t1, t2), (Δ1, Δ2))
[tr(h * σ) for σ in paulis]
h0 = substitute.(h, k => 0)
hpi = substitute.(h, k => pi)
[tr(h0 * σ) for σ in paulis]
[tr(hpi * σ) for σ in paulis]
##
function fast_energies(m::SMatrix{2,2})
    t = tr(m)
    d = det(m)
    s = sqrt(t^2 - 4d)
    (real(t + s) / 2, real(t - s) / 2)
end
function fast_energy_gap(m::SMatrix{2,2,Complex{Num},4})
    es = fast_energies(m)
    minimum(abs, es)
end
function LinearAlgebra.eigvals(m::SMatrix{2,2,Complex{Num},4})
    t = tr(m) |> simplify
    d = det(m) |> simplify
    [(t + sqrt(real(t^2 - 4d))) / 2 |> simplify,
        (t - sqrt(real(t^2 - 4d))) / 2 |> simplify]
end
bdgH(k, ε, (t * exp(1im * δϕ), 0), (Δ, 0)) |> eigvals
##
substitute.(skewH(k, ε, (t1, t2), (Δ1, Δ2); check=false), k => 0)
##

# function hoppingH(k, ε, (t1, t2))
#     t = 2 * real(cis(k) * t1 + cis(2k) * t2)
#     [ε+t 0; 0 ε+t]
# end
# function pairingH(k, (Δ1, Δ2))
#     Δ = 2 * 1im * imag(cis(k) * Δ1 + cis(2k) * Δ2)
#     [0 Δ; -Δ 0]
# end
#bdgH(k, ε, (t1, t2), (Δ1, Δ2); full=false, kwargs...) = (m = BdGMatrix(hoppingH(k, ε, (t1, t2)), pairingH(k, (Δ1, Δ2)); kwargs...); full ? m : m[2:3, 2:3])
skewH(k, ε, (t1, t2), (Δ1, Δ2); kwargs...) = QuantumDots.bdg_to_skew(bdgH(k, ε, (t1, t2), (Δ1, Δ2)); kwargs...)
function topoQ(ε, (t1, t2), (Δ1, Δ2); kwargs...)
    # pf1 = pfaffian(skewH(0, ε, (t1, t2), (Δ1, Δ2); kwargs...))
    # pf2 = pfaffian(skewH(pi, ε, (t1, t2), (Δ1, Δ2); kwargs...))
    # sign(pf1 * pf2)
    sign((2real(t1 + t2) + ε) * (2real(-t1 + t2) + ε))
end
sqrt(real(tr(bdgH(k, ε, (t1, t2), (Δ1, Δ2))^2)))

_energy_gap(k, p) = fast_energy_gap(bdgH(only(k), p...))
opt_prob = OptimizationProblem(_energy_gap, [0.0], ; lb=-pi, ub=pi)
function energy_gap(ε, (t1, t2), (Δ1, Δ2); nk=1000)
    # f(k, p) = abs(eigvals!(bdgH(only(k), ε, t, Δ), sortby=abs)[1])
    #f(k, p) = sqrt(real(tr(bdgH(only(k), ε, t, Δ)^2)))

    # f(_k, p) = (k = only(_k);
    # prob = OptimizationProblem(f, [0.0]; lb=-pi, ub=pi)
    # sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxtime)
    # f(k) = sqrt(real(tr(bdgH(k, ε, (t1, t2), (Δ1, Δ2))^2)))

    # function f(k)
    #     s1, c1 = sincos(k)
    #     s2, c2 = sincos(2k)
    #     # sqrt((-0.5ε - cos(2k)*real(t2) - cos(k)*real(t1) - sin(k)*imag(t1) - imag(t2)*sin(2k))^2 + ((1//2)*ε + cos(2k)*real(t2) + cos(k)*real(t1) - sin(k)*imag(t1) - imag(t2)*sin(2k))^2 + 2((-sin(k)*imag(Δ1) - imag(Δ2)*sin(2k))^2) - 2(-sin(k)*real(Δ1) - real(Δ2)*sin(2k))*(sin(k)*real(Δ1) + real(Δ2)*sin(2k)))

    #     sqrt((-0.5ε - c2 * real(t2) - c1 * real(t1) - s1 * imag(t1) - imag(t2) * s2)^2 + ((1 // 2) * ε + c2 * real(t2) + c1 * real(t1) - s1 * imag(t1) - imag(t2) * s2)^2 + 2((-s1 * imag(Δ1) - imag(Δ2) * s2)^2) - 2(-s1 * real(Δ1) - real(Δ2) * s2) * (s1 * real(Δ1) + real(Δ2) * s2))
    # end
    # f(k) = abs(det(bdgH(only(k), ε, (t1, t2), (Δ1, Δ2))))
    # f(k) = abs(eigvals(bdgH(only(k), ε, (t1, t2), (Δ1, Δ2)))[1])

    # p = (ε, (t1, t2), (Δ1, Δ2))
    # prob = remake(opt_prob; p)
    # sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited())
    # _energy_gap(sol, p)

    # f(k) = fast_energy_gap(bdgH(only(k), ε, (t1, t2), (Δ1, Δ2)))
    # ks = range(-pi, pi, nk)
    # Es = [f(k) for k in ks]
    # Emin, k = findmin(Es)
    # return Emin
    f1(k) = real(fast_energies(bdgH(k, ε, (t1, t2), (Δ1, Δ2)))[1])
    # f2(k) = real(fast_energies(bdgH(k, ε, (t1, t2), (Δ1, Δ2)))[2])
    Df1(x) = ForwardDiff.derivative(f1, x)
    # Df2(x) = ForwardDiff.derivative(f2, x)
    # fx = ZeroProblem(Df, 0)
    # sol = solve(Df, Roots.Newton())#, atol=0.0, rtol=0.0)
    sol1 = find_zeros(f1, -pi, pi)#, atol=0.0, rtol=0.0)
    if length(sol1) > 0
        return minimum(abs ∘ f1, sol1)
    end
    sol2 = find_zeros(Df1, -pi, pi)#, atol=0.0, rtol=0.0)
    return minimum(abs ∘ f1, sol2)

    # sol2 = find_zeros(f1, (0, pi))#, atol=0.0, rtol=0.0)
    # return minimum(f1, sol1)
    # return vcat(sol1, sol2)
    # return minimum(abs ∘ f1, vcat(sol1, sol2))
    # return map(f1, §vcat(sol1, sol2))
    # find_zero(f, 0, Order1())
    # ks = range(-pi, pi, nk)
    # Es = [f(k) for k in ks]
    # Emin, k = findmin(Es)
    # return Emin
end
##
@time energy_gap(0.5, (2exp(1im * pi / 4), 0), (1, 0))
@profview energy_gap(-1, (1 + 0.1im, 0), (1, 0))

##
let ks = range(-pi, pi, 100), f(k) = collect(fast_energies(bdgH(k, 0.5, (2exp(1im * pi / 4), 0), (1, 0)))), Df
    Df = x -> ForwardDiff.derivative(f, x)
    plot(ks, stack([Df(k) for k in ks])')
end |> display
##

h = bdgH(1, 0, (1 + 1im, 0), (1, 2))
skewH(pi, 0, (1 + 1im, 0), (1, 2); check=false)
[pfaffian(skewH(k, 0.2, (1, 0), (1, 0))) for k in (0, pi)]

##
@time [energy_gap(eps, (exp(1im * 1) * t, 0), (1, 0)) for eps in range(-4, 4, length=10), t in range(0, 4, length=11)];
@profview [energy_gap(eps, (exp(1im * 1) * t, 0), (1, 0)) for eps in range(-4, 4, length=10), t in range(0, 4, length=11)];

##
let ϕ = 0.5 * pi / 2
    ϵs = range(-4, 4, length=50)
    ts = range(0, 4, length=101)
    #ks = range(-pi, pi, length=100)
    data = [topoQ(eps, (exp(1im * ϕ) * t, 0), (1, 0); check=false) for eps in ϵs, t in ts]
    data2 = [energy_gap(eps, (exp(1im * ϕ) * t, 0), (1, 0)) for eps in ϵs, t in ts]
    p1 = heatmap(ts, ϵs, data, clims=(-1, 1), c=:redsblues, xlabel="t", ylabel="ε", title="Topological invariant")
    p2 = heatmap(ts, ϵs, data2, clims=(0, 1), c=:viridis, xlabel="t", ylabel="ε", title="Energy gap")
    plot(p1, p2) |> display
end
##
let t = 1, ϕ = 0.9 * pi / 2
    ϵs = range(-4, 4, length=100)
    Δs = range(-2, 2, length=101)
    #ks = range(-pi, pi, length=100)
    data = [topoQ(eps, (exp(1im * ϕ) * t, 0), (Δ, 0); check=false) for eps in ϵs, Δ in Δs]
    data2 = [energy_gap(eps, (exp(1im * ϕ) * t, 0), (Δ, 0)) for eps in ϵs, Δ in Δs]
    p1 = heatmap(ϵs, Δs, data', clims=(-1, 1), c=:redsblues, xlabel="ε", ylabel="Δ", title="Topological invariant")
    p2 = heatmap(ϵs, Δs, data2', clims=(0, 1), c=:viridis, xlabel="ε", ylabel="Δ", title="Determinant")
    plot(p1, p2) |> display
end
##
let t = 1
    ks = range(-pi, pi, length=100)
    energies = [[eigvals(bdgH(k, eps, (1, 0.0), (1.1, 0.0))) for k in ks] for eps in [-2t, 0, 2t]]
    plot([plot(ks, stack(es)') for es in energies]...)
end
##
let t = 1
    ks = range(-pi, pi, length=100)
    energies = [[eigvals(bdgH(k, eps, (1im, 0.0), (1.1, 0.0))) for k in ks] for eps in [-2t, 0, 2t]]
    plot([plot(ks, stack(es)') for es in energies]...)
end
##
bdgH(0.3, -1, (0.5, 0.0), (1, 0))
skewH(0, 0.2, (0.5, 0.0), (1, 0))

##
function effective_coeffs(; ε, t, θ, Δ, Ez, δϕ)
    Δ_aa, Δ_bb, Δ_ab, Δ_ba, t_aa, t_ab, t_ba, t_bb = LongerPoorMansMajoranas.perturbative_coeffs_homogeneous(; ε, t, θ, Δ, δϕ)
    ε0 = Ez - sqrt(Δ^2 + ε^2)
    (; E, ε_ba, ε_ab, t_nn, Δ_nn) = LongerPoorMansMajoranas.second_order_coeffs_homogeneous(; Δ_ba, Δ_ab, t_ba, t_ab, Ez, ε, Δ)
    (; ε=real(ε_ba + ε_ab) + ε0, t1=t_aa, t2=t_nn, Δ1=-conj(Δ_aa), Δ2=Δ_nn)
end
function topoQ(; first_order=false, ε, t, θ, Δ, Ez, δϕ)
    (; ε, t1, t2, Δ1, Δ2) = effective_coeffs(; ε, t, θ, Δ, Ez, δϕ)
    # (; ε, t1, t2, Δ1, Δ2) = effective_coeffs(; ε, t, θ, Δ, Ez, δϕ)
    # println("B:", ε)
    topoQ(ε, first_order ? (t1, 0t2) : (t1, t2), first_order ? (Δ1, 0Δ2) : (Δ1, Δ2); check=false)
end
function energy_gap(; first_order=false, ε, t, θ, Δ, Ez, δϕ, kwargs...)
    (; ε, t1, t2, Δ1, Δ2) = effective_coeffs(; ε, t, θ, Δ, Ez, δϕ)
    energy_gap(ε, first_order ? (t1, 0t2) : (t1, t2), first_order ? (Δ1, 0Δ2) : (Δ1, Δ2); kwargs...)
end
topoQ(; ε=1, δϕ=1, fixedparams...)
##
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), Δ=1, Ez=3)
ϵs = fixedparams.Ez - 0.1 .+ 0.4 .* range(-1, 1, length=150)
δϕs = range(0, pi, length=151)
##
dataE1 = [energy_gap(; ε, first_order=true, δϕ, fixedparams...) for ε in ϵs, δϕ in δϕs]
dataQ1 = [topoQ(; first_order=true, ε, δϕ, fixedparams...) for ε in ϵs, δϕ in δϕs]
dataE2 = [energy_gap(; ε, first_order=false, δϕ, fixedparams...) for ε in ϵs, δϕ in δϕs]
dataQ2 = [topoQ(; first_order=false, ε, δϕ, fixedparams...) for ε in ϵs, δϕ in δϕs]
##
p1E = heatmap(ϵs, δϕs, dataE1', clims=(1e-10, 1e-2), colorbar_scale=:identity, c=:viridis, xlabel="ε", ylabel="δϕ", title="Energy gap 1st order")
p1Q = heatmap(ϵs, δϕs, dataQ1', clims=(-1, 1), c=:redsblues, xlabel="ε", ylabel="δϕ", title="Topological invariant")
p2E = heatmap(ϵs, δϕs, dataE2', clims=(1e-10, 1e-2), colorbar_scale=:identity, c=:viridis, xlabel="ε", ylabel="δϕ", title="Energy gap 2nd order")
p2Q = heatmap(ϵs, δϕs, dataQ2', clims=(-1, 1), c=:redsblues, xlabel="ε", ylabel="δϕ", title="Topological invariant")
plot(p1E, p1Q, p2E, p2Q)
##
a = FermionBdGBasis(1:20)
c = FermionBdGBasis(1:20, (:↑, :↓))
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
p3 = heatmap(ϵs, δϕs, map(LDbdg, data_p1)', c=:viridis, xlabel="ε", ylabel="δϕ", title="LD, 2nd order, 20 sites")
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