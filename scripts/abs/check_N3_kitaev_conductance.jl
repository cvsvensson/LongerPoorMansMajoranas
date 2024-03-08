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
plot(μ1s, mapreduce(d -> first(d.energies), hcat, datal[:, 1, 1])', title="kitaev, μ12");
plot!(μ1s, mapreduce(d -> last(d.energies), hcat, datal[:, 1, 1])', ls=:dash)
plot(μ1s, mapreduce(d -> first(d.energies), hcat, datal3[:, 1, 1])', title="kitaev, μ123");
plot!(μ1s, mapreduce(d -> last(d.energies), hcat, datal3[:, 1, 1])', ls=:dash)
##
p = plot(; frame=:box, xlabel="ε", ylabel="δE", thickness_scaling=1.3, leg=:top)
for (data, l1, ls) in zip([csp..., cs], ["H1", "H2", "H"], [:solid, :dash, :dashdot])
    for (d, l2, c) in zip(data[[:onel, :twol, :threel]], ["ε1", "ε12", "ε123"], 1:3)
        plot!(p, εs, map(x -> (abs(x.gap)), d), label=string(l1, " ", l2); ls, c)
    end
end
p
##
pf_1 = plot(F3data["x"][1:6:end], mapreduce(d -> -first(d.energies)[1:4], hcat, cs.onel[:, 1])', title="full, μ1");
plot!(pf_1, F3data["x"][1:6:end], mapreduce(d -> -last(d.energies)[1:4], hcat, cs.onel[:, 1])', ls=:dash)
pf_12 = plot(F3data["x"][1:6:end], mapreduce(d -> -first(d.energies)[1:4], hcat, cs.twol[:, 1])', title="full, μ12");
plot!(pf_12, F3data["x"][1:6:end], mapreduce(d -> -last(d.energies)[1:4], hcat, cs.twol[:, 1])', ls=:dash)
pf_123 = plot(F3data["x"][1:6:end], mapreduce(d -> -first(d.energies)[1:4], hcat, cs.threel[:, 1])', title="full, μ123");
plot!(pf_123, F3data["x"][1:6:end], mapreduce(d -> -last(d.energies)[1:4], hcat, cs.threel[:, 1])', ls=:dash)
##
cspfdata = [csp[1:2]..., cs]
p1_1 = plot(F3data["x"][1:6:end], mapreduce(d -> first(d.energies)[1:4], hcat, csp[1].onel[:, 1])', title="p1, μ1");
plot!(p1_1, F3data["x"][1:6:end], mapreduce(d -> last(d.energies)[1:4], hcat, csp[1].onel[:, 1])', ls=:dash)
p1_12 = plot(F3data["x"][1:6:end], mapreduce(d -> first(d.energies)[1:4], hcat, csp[1].twol[:, 1])', title="p1, μ12");
plot!(p1_12, F3data["x"][1:6:end], mapreduce(d -> last(d.energies)[1:4], hcat, csp[1].twol[:, 1])', ls=:dash)
p1_123 = plot(F3data["x"][1:6:end], mapreduce(d -> first(d.energies)[1:4], hcat, csp[1].threel[:, 1])', title="p1, μ123");
plot!(p1_123, F3data["x"][1:6:end], mapreduce(d -> last(d.energies)[1:4], hcat, csp[1].threel[:, 1])', ls=:dash)
##
p2_1 = plot(F3data["x"][1:6:end], mapreduce(d -> first(d.energies)[1:4], hcat, csp[2].onel[:, 1])', title="p2, μ1");
plot!(p2_1, F3data["x"][1:6:end], mapreduce(d -> last(d.energies)[1:4], hcat, csp[2].onel[:, 1])', ls=:dash)
p2_12 = plot(F3data["x"][1:6:end], mapreduce(d -> first(d.energies)[1:4], hcat, csp[2].twol[:, 1])', title="p2, μ12");
plot!(p2_12, F3data["x"][1:6:end], mapreduce(d -> last(d.energies)[1:4], hcat, csp[2].twol[:, 1])', ls=:dash)
p2_123 = plot(F3data["x"][1:6:end], mapreduce(d -> first(d.energies)[1:4], hcat, csp[2].threel[:, 1])', title="p2, μ123");
plot!(p2_123, F3data["x"][1:6:end], mapreduce(d -> last(d.energies)[1:4], hcat, csp[2].threel[:, 1])', ls=:dash)
##
plot(pf_1, pf_12, pf_123, p1_1, p1_12, p1_123, p2_1, p2_12, p2_123, leg=false, size=(1200, 800))
##
Vs = [0.0]
T = 1 / 20
εs = range(2, 4, 100)
cs = conductance_sweep(F3data["fixedparams"], F3data["ss"].sol, εs, Vs, T)
csp = pert_conductance_sweep(F3data["fixedparams"], F3data["ss"].sol, εs, Vs, T)
##
pl = plot(μ1s, map(MPU, datal[:, 1, 1]), c=cgrad(:amp, rev=true), ylims=(0, 1), colorbar=false)
pl3 = heatmap(μ1s, Vls, map(MPU, datal3[:, :, 1]) |> permutedims, c=cgrad(:amp, rev=true), clims=(0, 1), colorbar=true)
pr = heatmap(μ1s, Vls, map(MPU, datar[:, :, 1]) |> permutedims, c=cgrad(:amp, rev=true), clims=(0, 1))
plot(pl, pl3)
