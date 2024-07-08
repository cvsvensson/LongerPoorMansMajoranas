using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra 
using CairoMakie
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
##
skewH(k, ε, (t1, t2), (Δ1, Δ2); kwargs...) = QuantumDots.bdg_to_skew(bdgH(k, ε, (t1, t2), (Δ1, Δ2)); kwargs...)
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
function effective_coeffs(; ε, t, θ, Δ, Ez, δϕ)
    Δ_aa, Δ_bb, Δ_ab, Δ_ba, t_aa, t_ab, t_ba, t_bb = LongerPoorMansMajoranas.perturbative_coeffs_homogeneous(; ε, t, θ, Δ, δϕ)
    ε0 = Ez - sqrt(Δ^2 + ε^2)
    (; E, ε_ba, ε_ab, t_nn, Δ_nn) = LongerPoorMansMajoranas.second_order_coeffs_homogeneous(; Δ_ba, Δ_ab, t_ba, t_ab, Ez, ε, Δ)
    (; ε=real(ε_ba + ε_ab) + ε0, t1=t_aa, t2=t_nn, Δ1=-conj(Δ_aa), Δ2=Δ_nn)
end
function topoQ(; first_order=false, ε, t, θ, Δ, Ez, δϕ, kwargs...)
    (; ε, t1, t2, Δ1, Δ2) = effective_coeffs(; ε, t, θ, Δ, Ez, δϕ)
    topoQ(ε, first_order ? (t1, 0t2) : (t1, t2), first_order ? (Δ1, 0Δ2) : (Δ1, Δ2); check=false)
end
function energy_gap(; first_order=false, ε, t, θ, Δ, Ez, δϕ, kwargs...)
    (; ε, t1, t2, Δ1, Δ2) = effective_coeffs(; ε, t, θ, Δ, Ez, δϕ)
    energy_gap(ε, first_order ? (t1, 0t2) : (t1, t2), first_order ? (Δ1, 0Δ2) : (Δ1, Δ2))
end
function bdgH(k; first_order=false, ε, t, θ, Δ, Ez, δϕ)
    (; ε, t1, t2, Δ1, Δ2) = effective_coeffs(; ε, t, θ, Δ, Ez, δϕ)
    bdgH(k, ε, first_order ? (t1, 0t2) : (t1, t2), first_order ? (Δ1, 0Δ2) : (Δ1, Δ2))
end
##
_data = load(datadir("final_data", "40-site-tuning.jld2"));
##
@unpack N, fixedparams, data, εs, δϕs = _data
##
dataE1 = [energy_gap(; ε, first_order=true, δϕ, fixedparams...) for ε in εs, δϕ in δϕs]
dataQ1 = [topoQ(; first_order=true, ε, δϕ, fixedparams...) for ε in εs, δϕ in δϕs]
@time dataE2 = [energy_gap(; ε, first_order=false, δϕ, fixedparams...) for ε in εs, δϕ in δϕs];
εs_high_res = range(first(εs), last(εs), length=500)
δϕs_high_res = range(first(δϕs), last(δϕs), length=500)
@time dataQ2 = [topoQ(; first_order=false, ε, δϕ, fixedparams...) for ε in εs_high_res, δϕ in δϕs_high_res];
##
cbwidth = 10
linewidth = 1.3
fig_phases = with_theme(theme_latexfonts()) do
    fig = Figure(size=0.8 .* (600, 300), fontsize=20, figure_padding=5)

    g = fig[1, 1] = GridLayout()
    xticks = WilkinsonTicks(3)
    ax = Axis(g[1, 1]; xlabel=L"ε/Δ", title=L"N=%$N", ylabel=L"δϕ", yticks=(pi * [0, 1 / 2, 1], [L"0", L"\frac{\pi}{2}", L"π"]), xticks)
    ax1 = Axis(g[1, 3]; xlabel=L"ε/Δ", title=L"N = \infty, \,\, H_2^\text{eff}", ylabel=L"δϕ", yticks=(pi * [0, 1 / 2, 1], [L"0", L"\frac{\pi}{2}", L"π"]), xticks)

    linkaxes!(ax, ax1)
    hideydecorations!(ax1)

    target = x -> LD_cells(x)

    hmap_kwargs = (; colormap=Reverse(:viridis), colorscale=identity, colorrange=(0, maximum(map(target, data))))
    contour_kwargs = (; color=:red, levels=[0.0], linewidth=1.3)

    hmap = heatmap!(ax, εs, δϕs, map(target, data)'; hmap_kwargs...)
    f_egap = heatmap!(ax1, εs, δϕs, dataE2; colormap=Reverse(:davos), colorscale=log10, colorrange=(1e-3, 1e-1))
    f_Q = contour!(ax1, εs_high_res, δϕs_high_res, dataQ2; contour_kwargs...)
    l_Q = lines!(ax1, Float64[], Float64[]; label="Q=0", contour_kwargs...)
    axislegend(ax1; position=:lt)

    ticks = ([0, 0.25, 1 / 2, 1], ["0", ".25", "0.5", "1"])
    ticklabelsize = 16
    Colorbar(g[1, 2], hmap; width=cbwidth, ticksize=cbwidth, tickalign=true, ticks=LinearTicks(3), ticklabelsize)
    Label(g[1, 2, Bottom()], " LD", tellwidth=false, tellheight=false, fontsize=20)
    Label(g[1, 4, Bottom()], "   |δE/Δ|", tellwidth=false, tellheight=false, fontsize=20)

    Colorbar(g[1, 4], f_egap; width=cbwidth, ticksize=cbwidth, tickalign=true, ticks=LogTicks(WilkinsonTicks(3)), ticklabelsize)
    colgap!(g, 1, 10)
    colgap!(g, 3, 10)
    fig
end
##
save(plotsdir("phase_diagrams_$(N)_vs_inf.png"), fig_phases; px_per_unit=10)
