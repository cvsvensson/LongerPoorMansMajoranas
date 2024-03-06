using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra#, BlackBoxOptim
using Plots
using Symbolics
using Folds
using Accessors
using JLD2
using Statistics
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
    H = blockdiagonal(0.0(c[1]' * c[3] + c[1]' * c[3]' + hc) + QuantumDots.kitaev_hamiltonian(c; μ=μ, fixedparams...), c)
    sol = fullsolve(H, c; transport)
    return sol
end

##
μ1s = range(-10, 10, 52)
μ2s = μ1s
μ3 = 0.0
Vls = range(-10, 10, 52)
ϕs = let n = 15
    range(0, 2π, n + 1)[1:n]
end
fixedparams = (; t=[1.0, 2.0, 0], V=0)
T = 1 / 10
# datal = [get_sol((μ1, μ1, μ1), (@insert fixedparams.Δ = fixedparams.t .* [1, exp(1im * ϕ), 0]), Transport(QuantumDots.PauliSystem, (; T, μ=(Vl, 0.0)))) for μ1 in μ1s, Vl in Vls, ϕ in ϕs]
datal = [get_sol((μ1, μ1, μ3), (@insert fixedparams.Δ = fixedparams.t .* [1, exp(1im * ϕ), 0]), Transport(QuantumDots.PauliSystem, (; T, μ=(Vl, 0.0)))) for μ1 in μ1s, Vl in Vls, ϕ in ϕs]
datal3 = [get_sol((μ1, μ1, μ1 + μ3), (@insert fixedparams.Δ = fixedparams.t .* [1, exp(1im * ϕ), 0]), Transport(QuantumDots.PauliSystem, (; T, μ=(Vl, 0)))) for μ1 in μ1s, Vl in Vls, ϕ in ϕs]
#datar = [get_sol((μ1, μ1, μ3), (@insert fixedparams.Δ = fixedparams.t .* [1, exp(1im * ϕ), 0]), Transport(QuantumDots.PauliSystem, (; T, μ=(0.0, Vl)))) for μ1 in μ1s, Vl in Vls, ϕ in ϕs]
# datar = [get_sol((μ12..., μ3), fixedparams, Transport(QuantumDots.PauliSystem, (; T=1 / 10, μ=(0.0, Vl)))) for μ12 in zip(μ1s, μ2s), Vl in Vls]
##
cl = mean(map(d -> d.conductance[1, 1], datal), dims=3)[:, :, 1]
cl3 = mean(map(d -> d.conductance[1, 1], datal3), dims=3)[:, :, 1]
cr = mean(map(d -> d.conductance[2, 2], datar), dims=3)[:, :, 1]
##
colorscale = cgrad(:amp, scale=:exp)#cgrad(:seaborn_rocket_gradient, rev=true)
pl = heatmap(μ1s, Vls, cl |> permutedims, c=colorscale, clims=(0, 5), colorbar=false)
pl3 = heatmap(μ1s, Vls, cl3 |> permutedims, c=colorscale, clims=(0, 5), colorbar=false)
#pr = heatmap(μ1s, Vls, cr |> permutedims, c=colorscale, clims=(0, 5))
plot(pl, pl3)
##
plot(μ1s, mapreduce(d -> first(d.energies), hcat, datal[:, 1, 1])', title = "kitaev, μ12");
plot!(μ1s, mapreduce(d -> last(d.energies), hcat, datal[:, 1, 1])', ls=:dash)
plot(μ1s, mapreduce(d -> first(d.energies), hcat, datal3[:, 1, 1])', title = "kitaev, μ123");
plot!(μ1s, mapreduce(d -> last(d.energies), hcat, datal3[:, 1, 1])', ls=:dash)
##
plot(F3data["x"][1:6:end], mapreduce(d -> -first(d.energies)[1:4], hcat, cs.twol[:, 1])', title = "full, μ12");
plot!(F3data["x"][1:6:end], mapreduce(d -> -last(d.energies)[1:4], hcat, cs.twol[:, 1])', ls=:dash)
plot(F3data["x"][1:6:end], mapreduce(d -> -first(d.energies)[1:4], hcat, cs.threel[:, 1])', title = "full, μ123");
plot!(F3data["x"][1:6:end], mapreduce(d -> -last(d.energies)[1:4], hcat, cs.threel[:, 1])', ls=:dash)

##
Vs = [0.0]
cs = conductance_sweep(F3data["fixedparams"], F3data["ss"].sol, F3data["x"][1:6:end], Vs, T)
##
pl = plot(μ1s, map(MPU, datal[:, 1, 1]), c=cgrad(:amp, rev=true), ylims=(0, 1), colorbar=false)
pl3 = heatmap(μ1s, Vls, map(MPU, datal3[:, :, 1]) |> permutedims, c=cgrad(:amp, rev=true), clims=(0, 1), colorbar=true)
pr = heatmap(μ1s, Vls, map(MPU, datar[:, :, 1]) |> permutedims, c=cgrad(:amp, rev=true), clims=(0, 1))
plot(pl, pl3)
