using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra#, BlackBoxOptim
using Plots
using Symbolics
using Folds
using Accessors
using JLD2
includet(scriptsdir("abs", "phase_plots.jl"))
includet(scriptsdir("abs", "phase_misc.jl"))

##
c = FermionBasis(1:3, qn=QuantumDots.parity)
transport = Transport(QuantumDots.PauliSystem, (; T=1 / 40, μ=(0.0, 0.0)))
##
fixedparams = (; t=1.0, Δ=1.0, V=0)
H = blockdiagonal(QuantumDots.kitaev_hamiltonian(c; μ=0, fixedparams...), c)
sol = fullsolve(H, c; transport)
##
function get_sol(μ, fixedparams, transport)
    H = blockdiagonal(0.2(c[1]' * c[3] + c[1]' * c[3]' + hc) + QuantumDots.kitaev_hamiltonian(c; μ=μ, fixedparams...), c)
    sol = fullsolve(H, c; transport)
    return sol
end

##
μ1s = range(-3, 3, 22)
μ2s = μ1s
μ3 = 0
Vls = range(-3, 3, 42)
ϕs = let n = 5
    range(0, 2π, n + 1)[1:n]
end
fixedparams = (; t=[1.0, 2.0, 0], V=0)
T = 1 / 20
# datal = [get_sol((μ1, μ1, μ1), (@insert fixedparams.Δ = fixedparams.t .* [1, exp(1im * ϕ), 0]), Transport(QuantumDots.PauliSystem, (; T, μ=(Vl, 0.0)))) for μ1 in μ1s, Vl in Vls, ϕ in ϕs]
datal = [get_sol((μ1, μ1, μ3), (@insert fixedparams.Δ = fixedparams.t .* [1, exp(1im * ϕ), 0]), Transport(QuantumDots.PauliSystem, (; T, μ=(Vl, 0.0)))) for μ1 in μ1s, Vl in Vls, ϕ in ϕs]
datal3 = [get_sol((μ1, μ1, μ1), (@insert fixedparams.Δ = fixedparams.t .* [1, exp(1im * ϕ), 0]), Transport(QuantumDots.PauliSystem, (; T, μ=(Vl, 0)))) for μ1 in μ1s, Vl in Vls, ϕ in ϕs]
#datar = [get_sol((μ1, μ1, μ3), (@insert fixedparams.Δ = fixedparams.t .* [1, exp(1im * ϕ), 0]), Transport(QuantumDots.PauliSystem, (; T, μ=(0.0, Vl)))) for μ1 in μ1s, Vl in Vls, ϕ in ϕs]
# datar = [get_sol((μ12..., μ3), fixedparams, Transport(QuantumDots.PauliSystem, (; T=1 / 10, μ=(0.0, Vl)))) for μ12 in zip(μ1s, μ2s), Vl in Vls]
##
using Statistics
cl = mean(map(d -> d.conductance[1, 1], datal), dims=3)[:, :, 1]
cl3 = mean(map(d -> d.conductance[1, 1], datal3), dims=3)[:, :, 1]
cr = mean(map(d -> d.conductance[2, 2], datar), dims=3)[:, :, 1]
##
colorscale = :amp#cgrad(:seaborn_rocket_gradient, rev=true)
pl = heatmap(μ1s, Vls, cl |> permutedims, c=colorscale, clims=(0, 5), colorbar=false)
pl3 = heatmap(μ1s, Vls, cl3 |> permutedims, c=colorscale, clims=(0, 5), colorbar=false)
pr = heatmap(μ1s, Vls, cr |> permutedims, c=colorscale, clims=(0, 5))
plot(pl, pl3)

##
pl = plot(μ1s, map(MPU, datal[:, 1, 1]), c=cgrad(:amp, rev=true), ylims=(0, 1), colorbar=false)
pl3 = heatmap(μ1s, Vls, map(MPU, datal3[:, :, 1]) |> permutedims, c=cgrad(:amp, rev=true), clims=(0, 1), colorbar=true)
pr = heatmap(μ1s, Vls, map(MPU, datar[:, :, 1]) |> permutedims, c=cgrad(:amp, rev=true), clims=(0, 1))
plot(pl, pl3)
