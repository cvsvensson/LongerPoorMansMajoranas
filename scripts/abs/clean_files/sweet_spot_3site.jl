using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra#, BlackBoxOptim
using Plots
using Symbolics
using Folds
using Accessors
using JLD2
using DataFrames
using LaTeXStrings
using FiniteDiff
synceddir(args...) = joinpath(ENV["Dropbox"], "data", "LongerPoorMans", args...)
includet(scriptsdir("abs", "phase_plots.jl"))
includet(scriptsdir("abs", "phase_misc.jl"))
##
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0, Ez=3)
## 2site sweet spot
bdg = iszero(fixedparams.V) && iszero(fixedparams.U)
target = bdg ? LDbdg : LD
c2 = bdg ? FermionBdGBasis(1:2, (:↑, :↓)) : FermionBasis(1:2, (:↑, :↓); qn=QuantumDots.parity)
f2, f2!, cache2 = hamfunc(Hδϕ_Hε(), c2, fixedparams)
exps = range(0.1, 4, 5)
_homogeneous_ss2 = find_sweet_spot((f2, f2!, cache2), c2, Hδϕ_Hε(); exps, MaxTime=1, target, minexcgap=0, alg=BestOf(best_algs()[1:end-1]))
homogeneous_ss2 = merge(_homogeneous_ss2, get_gap_derivatives(f2, c2, _homogeneous_ss2.sol))

##
initials = [2.2, sqrt(fixedparams.Ez^2 - fixedparams.Δ^2)]#2.95]
N = 3
c = bdg ? FermionBdGBasis(1:N, (:↑, :↓)) : FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
f, f!, cache = hamfunc(Rδϕ_Rε(), c, fixedparams)
fh, fh!, cacheh = hamfunc(Hδϕ_Hε(), c, fixedparams)

f([1, 1, 1]) - fh([1, 1]) |> norm

exps = range(0.1, 4, 5)

param_cost = (ps, p, exp) -> 10.0^(-exp) * 1 / (10.0^(-exp) + sum(p0 -> norm(p0 .- p)^4, ps))
alg = LongerPoorMansMajoranas.MultipleSS(BestOf(best_algs()[1:end-1]), 2, param_cost)
# alg = (BestOf(best_algs()[1:end-1]))
_homogeneous_ss = find_sweet_spot((fh, fh!, cacheh), c, Hδϕ_Hε(); exps, MaxTime=1, target, minexcgap=0, alg, initials)
homogeneous_ss = map(sol -> merge(sol, get_gap_derivatives(fh, c, sol.sol)), _homogeneous_ss)
# δϕ0 = homogeneous_ss.sol[1]

##
ss_sorter(sols) = sort(sols, by=x -> x.sol[1])
##
ss = ss_sorter(homogeneous_ss)[1]
ε0 = ss.sol[2]
initials = collect(ss.sol)
MaxTime = 1
minexcgap = ss.optsol.excgap * 0.0
phase_data = []
level_data = []
δε2s = range(-0.2, 0.4, length=20) #Phase branch
# δε2s = range(-0.2, 0.7, length=10) #level branch
for δε2 in δε2s
    # alg = BestOf(best_algs()[1:end-1])
    # alg = LongerPoorMansMajoranas.MultipleSS(BestOf(best_algs()[1:end-1]), 2, param_cost)
    alg = BestOf(best_algs()[1:end-1])
    hamfunc = δϕε1 -> f!(cache, [δϕε1[1], δϕε1[2], δϕε1[2] + δε2])
    prob = OptProb(; hamfunc, basis=c, optparams=Rδϕ_Rε(), target)

    #WARNING: This is a hack to get the phase branch
    phase_ranges = [(0.0, 2.0), ε0 .+ 1 .* (-0.5, 0.05)] #phase branch
    level_ranges = [(0.0, 1.0pi), ε0 .+ 1 .* (0.025, 0.3)] # level branch

    # ranges = [(0.0, 1.0pi), ε0 .+ 2 .* (-1, 1)]

    # println("INITIALS: ", initials)
    # sols = solve(prob, alg; minexcgap, maxiters=10000, MaxTime, exps, initials, ranges)
    # phase branch
    phase_sols = [solve(prob, alg; minexcgap, maxiters=10000, MaxTime, exps, initials, ranges=phase_ranges)]
    phase_sols = ss_sorter(phase_sols)
    phase_sols2 = map(sol -> merge(sol, get_gap_derivatives(hamfunc, c, sol.sol)), phase_sols)
    push!(phase_data, phase_sols2)
    # repeat for level branch
    level_sols = [solve(prob, alg; minexcgap, maxiters=10000, MaxTime, exps, initials, ranges=level_ranges)]
    level_sols = ss_sorter(level_sols)
    level_sols2 = map(sol -> merge(sol, get_gap_derivatives(hamfunc, c, sol.sol)), level_sols)
    push!(level_data, level_sols2)
end
##
let phase_data = map(x -> x[1], phase_data), level_data = map(x -> x[1], level_data)
    p_LD = plot(δε2s, map(x -> LDbdg(x.optsol), phase_data), markers=true, ylabel="LD", xlabel="δε2", label="Inhomogeneous phase branch", ylims=(0, 1))
    plot!(δε2s, map(x -> LDbdg(x.optsol), level_data), markers=true, ylabel="LD", xlabel="δε2", label="Inhomogeneous level branch")
    # scatter!([0], [LDbdg(ss.optsol)], c=:red, label="homogeneous")
    hline!([LDbdg(homogeneous_ss2.optsol)], label="2 site sweet spot", lw=2, ls=:dash)
    p_excgap = plot(δε2s, map(x -> x.optsol.excgap, phase_data), markers=true, ylabel="excgap", xlabel="δε2", label="Inhomogeneous phase branch", ylims=(0, 0.3))
    plot!(δε2s, map(x -> x.optsol.excgap, level_data), markers=true, ylabel="excgap", xlabel="δε2", label="Inhomogeneous level branch")
    # scatter!([0], [ss.optsol.excgap], c=:red, label="homogeneous")
    hline!([homogeneous_ss2.optsol.excgap], ls=:dash, lw=2, label="2 site sweet spot")
    p_gap = plot(δε2s, map(x -> x.optsol.gap, phase_data), markers=true, ylabel="gap", xlabel="δε2", label="Inhomogeneous phase branch")
    plot!(δε2s, map(x -> x.optsol.gap, level_data), markers=true, ylabel="gap", xlabel="δε2", label="Inhomogeneous level branch")
    # scatter!([0], [ss.optsol.gap], c=:red, label="homogeneous")
    hline!([homogeneous_ss2.optsol.gap], ls=:dash, lw=2, label="2 site sweet spot")
    p_gap_der = plot(δε2s, map(x -> norm(x.gradient), phase_data), markers=true, ylabel="gap_der", xlabel="δε2", label="Inhomogeneous phase branch", ylims=(0, 1))
    plot!(δε2s, map(x -> norm(x.gradient), level_data), markers=true, ylabel="gap_der", xlabel="δε2", label="Inhomogeneous level branch", ylims=(0, 1))
    # scatter!([0], [norm(ss.gradient)], c=:red, label="homogeneous")
    hline!([norm(homogeneous_ss2.gradient)], ls=:dash, lw=2, label="2 site sweet spot")

    p_params = plot(δε2s, stack(map(x -> (x.sol), phase_data))', markers=true, ylabel="params", xlabel="δε2", ylims=(0, pi))
    plot!(δε2s, stack(map(x -> (x.sol), level_data))', markers=true, ylabel="params", xlabel="δε2", ylims=(0, pi))
    # scatter!([0], [norm(ss.gradient)], c=:red, label="homogeneous")

    # plot(p_LD, p_excgap, p_params)#, layout=(1, 2))
    plot(p_LD, p_excgap, layout=(1, 2))
    # plot(p_LD, p_excgap, p_gap_der, p_gap, p_params, layout=(3, 2), size=(800, 800))
end
##
savefig("inhomogeneous_sweet_spot_3site.png")
##
for n in [3, 5, 10] # Look at tuning plot for these detunings
    let k = n, hamfunc, a, data_p1, data = map(x -> x[1], data)
        a = FermionBdGBasis(1:3)
        δε2 = δε2s[k]
        δϕ, ε1 = data[k].sol
        hamfunc = (δϕ, ε1) -> f!(cache, [δϕ, ε1, ε1 + δε2])
        εs = ε1 .+ 0.2 .* range(-1, 1, 50)
        δϕs = range(0, pi, 51)
        newdata = [(
            begin
                H = hamfunc(δϕ, ε1)
                merge(fullsolve(H, c),
                    get_gap_derivatives(x -> hamfunc(x...), c, [δϕ, ε1]))
            end
        ) for ε1 in εs, δϕ in δϕs] |> permutedims
        d = Dict("ss" => data[k], "data" => newdata, "N" => N, "labels" => ("ε1", "δϕ"), "x" => εs, "y" => δϕs)
        p1_hamfunc = (δϕ, ε1) -> LongerPoorMansMajoranas.perturbative_hamiltonian(a, 1; ε=[ε1, ε1 + δε2, ε1], δϕ=δϕ * [1, 1], t=fixedparams.t, θ=fixedparams.θ, Δ=fixedparams.Δ * [1, 1, 1], Ez=fixedparams.Ez)
        data_p1 = [(
            begin
                #H = LongerPoorMansMajoranas.perturbative_hamiltonian(a, 1; ε=[ε1, ε1 + δε2, ε1], δϕ=δϕ * [1, 1], t=fixedparams.t, θ=fixedparams.θ, Δ=fixedparams.Δ * [1, 1, 1], Ez=fixedparams.Ez)
                # H = LongerPoorMansMajoranas.perturbative_hamiltonian_homogeneous(a, 2; ε, δϕ, fixedparams...)
                merge(fullsolve(p1_hamfunc(δϕ, ε1), a),
                    get_gap_derivatives(x -> p1_hamfunc(x...), a, [δϕ, ε1]))
            end
        ) for ε1 in εs, δϕ in δϕs] |> permutedims

        data_coeffs = [(
            let
                Δ = fixedparams.Δ * [1, 1, 1]
                ε = [ε1, ε1 + δε2, ε1]
                θ = fixedparams.θ
                δϕ = δϕ * [1, 1]
                t = fixedparams.t
                Ez = fixedparams.Ez
                coeffs = stack(LongerPoorMansMajoranas.perturbative_coeffs(n; Δ, ε, θ, δϕ, t) for n in 1:2)
                Δ_aa, Δ_bb, Δ_ab, Δ_ba, t_aa, t_ab, t_ba, t_bb = collect(eachrow(coeffs))
                E1, E2, ε_ba, ε_ab, t_nn, Δ_nn = LongerPoorMansMajoranas.second_order_coeffs(1; Δ_ba, Δ_ab, t_ba, t_ab, Ez, ε, Δ)
                ε0 = [real(Ez - sqrt(Δ[n]^2 + ε[n]^2)) for n in 1:3]
                (; Δ_aa, Δ_bb, Δ_ab, Δ_ba, t_aa, t_ab, t_ba, t_bb, E1, E2, ε_ba, ε_ab, t_nn, Δ_nn, ε0)
            end
        ) for ε1 in εs, δϕ in δϕs] |> permutedims

        dp1 = Dict("ss" => data[k], "data" => data_p1, "N" => 3, "labels" => ("ε1", "δϕ"), "x" => εs, "y" => δϕs)
        p1 = plot_f(d, x -> LDbdg(x), clim_max=1, c=cgrad(:viridis, rev=true), ss_label="3-site sweet spot", legend=false)
        p2 = plot_f(dp1, x -> LDbdg(x), clim_max=1, c=cgrad(:viridis, rev=true), ss_label="3-site sweet spot", legend=:bottomright)
        # p1 = plot_f(d, x -> norm(x.gradient), clim_max=1, c=cgrad(:viridis, rev=true), ss_label="3-site sweet spot", legend=false)
        # p2 = plot_f(dp1, x -> norm(x.gradient), clim_max=1, c=cgrad(:viridis, rev=true), ss_label="3-site sweet spot", legend=:bottomright)
        # p1 = plot_f(d, x -> norm(x.gradient[2]) / LDf(x), clim_max=2, c=:redsblues, ss_label="3-site sweet spot", legend=false)
        # p2 = plot_f(dp1, x -> norm(x.gradient[2]) / LDf(x), clim_max=2, c=:redsblues, ss_label="3-site sweet spot", legend=:bottomright)
        # p1 = plot_f(d, x -> x.excgap, clim_max=.25, c=:viridis, ss_label="3-site sweet spot", legend = false)
        # p2 = plot_f(dp1, x -> x.excgap, clim_max=.25, c=:viridis, ss_label="3-site sweet spot")
        p3 = heatmap(εs, δϕs, map(x -> abs(x.t_aa[1]) - abs(x.Δ_aa[1]), data_coeffs), c=:redsblues, clims=0.1 .* (-1, 1))
        p4 = heatmap(εs, δϕs, map(x -> abs(x.Δ_nn) * angle(x.Δ_nn), data_coeffs), c=:redsblues, clims=0.001 .* (-1, 1))
        p5 = heatmap(εs, δϕs, map(x -> abs(x.t_nn) - abs(x.Δ_nn), data_coeffs), c=:redsblues, clims=0.1 .* (-1, 1))
        p6 = heatmap(εs, δϕs, map(x -> (x.ε_ba + x.ε_ab + x.ε0[2])^2 + (x.ε_ab + x.ε0[1])^2, data_coeffs), c=:redsblues, clims=0.01 .* (0, 1))


        # p = heatmap(εs, δϕs, map(x -> LDbdg(x), newdata)', c=:viridis, clims=(0, 1), xlabel="ε1", ylabel="δϕ", title="LD, 3 site, inhomogeneous")

        scatter!(p3, [ε1], [δϕ], c=:red)
        scatter!(p4, [ε1], [δϕ], c=:red)
        scatter!(p5, [ε1], [δϕ], c=:red)
        scatter!(p6, [ε1], [δϕ], c=:red)
        # p
        plot(p1, p2, p3, p4, p5, p6, plot_title="δε2 = $δε2", layout=(3, 2), size=(800, 1200)) |> display
    end
end
##

##
sweet_spots_Ez3 = []
sweet_spots_Ez2 = []
Ezs = range(1.1, 10, 10)
for Ez in Ezs
    fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0, Ez)
    f2, f2!, cache2 = hamfunc(Hδϕ_Hε(), c2, fixedparams)
    exps = range(0.1, 3, 4)
    MaxTime = 0.5
    initials = [2.2, sqrt(fixedparams.Ez^2 - fixedparams.Δ^2)]#2.95]
    homogeneous_ss2 = find_sweet_spot((f2, f2!, cache2), c2, Hδϕ_Hε(); exps, MaxTime, target, minexcgap=0, alg=BestOf(best_algs()[1:end-1]), initials)
    f, f!, cache = hamfunc(Rδϕ_Rε(), c, fixedparams)
    fh, fh!, cacheh = hamfunc(Hδϕ_Hε(), c, fixedparams)
    homogeneous_ss = find_sweet_spot((fh, fh!, cacheh), c, Hδϕ_Hε(); exps, MaxTime, target, minexcgap=0, alg=BestOf(best_algs()[1:end-1]), initials)
    # δϕ0 = homogeneous_ss.sol[1]
    push!(sweet_spots_Ez3, merge(homogeneous_ss, get_gap_derivatives(fh, c, homogeneous_ss.sol)))
    push!(sweet_spots_Ez2, merge(homogeneous_ss2, get_gap_derivatives(f2, c2, homogeneous_ss2.sol)))
end
##
p_Ez_excgap = plot(Ezs, [x.optsol.excgap for x in sweet_spots_Ez3], markers=true, ylabel="excgap", xlabel="Ez", ylims=(0, 0.8))
plot!(Ezs, [LDbdg(x.optsol) for x in sweet_spots_Ez2], markers=true, ls=:dash)
p_Ez_LD = plot(Ezs, [LDbdg(x.optsol) for x in sweet_spots_Ez3], markers=true, ylabel="LDbdg", xlabel="Ez", clims=(0, 1.5))
plot!(Ezs, [LDbdg(x.optsol) for x in sweet_spots_Ez2], markers=true, ls=:dash)
plot(p_Ez_excgap, p_Ez_LD)


##
inhomo_sols = []
for k in eachindex(Ezs)
    fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0, Ez=Ezs[k])
    alg = BestOf(best_algs()[1:end-1])
    f, f!, cache = hamfunc(Rδϕ_Rε(), c, fixedparams)
    minexcgap = sweet_spots_Ez3[k].optsol.excgap * 0.1
    δϕ0, ε0 = sweet_spots_Ez3[k].sol |> collect
    prob = OptProb(; hamfunc=x -> f!(cache, [x[1], x[2], x[2] + x[3]]), basis=c, optparams=Rδϕ_Rε(), target)
    ranges = [(0.0, 1.0pi), ε0 .+ 1 .* (-1, 1), (-0.5, 0.5)]
    initials = [δϕ0, ε0, 0.0]
    sol = solve(prob, alg; minexcgap, maxiters=10000, MaxTime=0.5, exps, initials, ranges)
    push!(inhomo_sols, merge(sol, get_gap_derivatives(f, c, sol.sol)))
end
##
c_main = 1
p_Ez_excgap = plot(Ezs, [x.optsol.excgap for x in sweet_spots_Ez3], markers=true, ylabel="excgap", xlabel="Ez", ylims=(0, 0.8), c=c_main, label="3-site homogeneous")
plot!(Ezs, [LDbdg(x.optsol) for x in sweet_spots_Ez2], markers=true, label="2-site")
plot!(Ezs, [LDbdg(x.optsol) for x in inhomo_sols], markers=true, ls=:dash, c=c_main, label="3-site inhomogeneous")
p_Ez_LD = plot(Ezs, [LDbdg(x.optsol) for x in sweet_spots_Ez3], markers=true, ylabel="LDbdg", xlabel="Ez", clims=(0, 1.5), c=c_main, label="3-site homogeneous")
plot!(Ezs, [LDbdg(x.optsol) for x in sweet_spots_Ez2], markers=true, label="2-site")
plot!(Ezs, [LDbdg(x.optsol) for x in inhomo_sols], markers=true, ls=:dash, c=c_main, label="3-site inhomogeneous")
p_Ez_gap = plot(Ezs, [abs(x.optsol.gap) for x in sweet_spots_Ez3], markers=true, ylabel="gap", xlabel="Ez", clims=(1e-5, 0.1), yscale=:log, c=c_main)
plot!(Ezs, [abs(x.optsol.gap) for x in sweet_spots_Ez2], markers=true)
plot!(Ezs, [abs(x.optsol.gap) for x in inhomo_sols], markers=true, ls=:dash, c=c_main)
p_Ez_grad = plot(Ezs, [norm(x.gradient) for x in sweet_spots_Ez3], markers=true, ylabel="|gradient|", xlabel="Ez", clims=(0, 1.5), c=c_main)
plot!(Ezs, [norm(x.gradient) for x in sweet_spots_Ez2], markers=true)
plot!(Ezs, [norm(x.gradient) for x in inhomo_sols], markers=true, ls=:dash, c=c_main)
p_Ez_hess = plot(Ezs, [norm(x.hessian) for x in sweet_spots_Ez3], markers=true, ylabel="|hessian|", xlabel="Ez", clims=(0, 1.5), c=c_main)
plot!(Ezs, [norm(x.hessian) for x in sweet_spots_Ez2], markers=true)
plot!(Ezs, [norm(x.hessian) for x in inhomo_sols], markers=true, ls=:dash, c=c_main)
plot(p_Ez_excgap, p_Ez_LD, p_Ez_grad, p_Ez_hess, p_Ez_gap)
plot(p_Ez_LD, p_Ez_excgap)
##
savefig("inhomogeneous_Vz_improvement.png")

##
plot(Ezs, stack(map((h, i) -> [h.sol[1] - i.sol[1], h.sol[2] - i.sol[2], h.sol[2] - i.sol[3] - i.sol[2]], sweet_spots_Ez3, inhomo_sols))', marker=true, label=["δϕ" "ε1" "ε2"], ylabel="difference", xlabel="Ez")
##
[LDbdg(fullsolve(hamfunc(Rδϕ_Rε(), c, merge(fixedparams, (; Ez)))[1](collect(sol.sol)), c)) for (Ez, sol) in zip(Ezs, inhomo_sols)] |> plot
[LDbdg(fullsolve(hamfunc(Hδϕ_Hε(), c, merge(fixedparams, (; Ez)))[1](collect(sol.sol)[1:2]), c)) for (Ez, sol) in zip(Ezs, inhomo_sols)] |> plot!
[LDbdg(fullsolve(hamfunc(Rδϕ_Rε(), c, merge(fixedparams, (; Ez)))[1]([sol.sol[1:2]..., sol.sol[3]]), c)) for (Ez, sol) in zip(Ezs, inhomo_sols)] |> plot!


##
let sol = inhomo_sols[1]
    get_gap_derivatives(hamfunc(Rδϕ_Rε(), c, merge(fixedparams, (; Ez=3)))[1], c, [sol.sol[1:2]..., sol.sol[3]])
end