using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra#, BlackBoxOptim
using Plots
using Symbolics
using Folds
using Accessors
using JLD2
includet(scriptsdir("abs", "phase_plots.jl"))
kitaev_ham(c, ε, Δ, t) = BdGMatrix(QuantumDots.kitaev_hamiltonian(c; μ=-ε, t, Δ); check=false)
##
##
N = 3
cK2 = FermionBdGBasis(1:2)
cKN = FermionBdGBasis(1:N)
kitaev_ham(c, ε, Δ, t) = BdGMatrix(QuantumDots.kitaev_hamiltonian(c; μ=-ε, t, Δ); check=false)
@variables t ε
kitaev_ham_fast2, kitaev_ham_fast2! = build_function(kitaev_ham(cK2, ε, 1, t), [t, ε], expression=Val{false})
kitaev_ham_fastN, kitaev_ham_fastN! = build_function(kitaev_ham(cKN, ε, 1, t), [t, ε], expression=Val{false})
cache2 = kitaev_ham_fast2([0.1, 0.1])
cacheN = kitaev_ham_fastN([0.1, 0.1])
##
@time h1 = kitaev_ham(cKN, 0.1, 1, 0.1);
@time h2 = kitaev_ham_fastN([0.1, 0.1]);
@time h3 = kitaev_ham_fastN!(cacheN, [0.1, 0.1]);
@assert h1 ≈ h2 ≈ h3
##
ts = range(-2, 2, length=100)
Kεs = range(-3, 3, length=100)
iter = Iterators.product(ts, Kεs) |> collect
# @time data = map(tε -> fullsolve(kitaev_ham(c, tε[2], 1, tε[1]), c), iter);
#@time data2 = map(tε -> fullsolve(kitaev_ham_fast2(tε), cK2), iter);
@time dataN = map(tε -> fullsolve(kitaev_ham_fastN!(cacheN, tε), cKN), iter);
##
kitaev_ss = (-1, 0)
kitaev_phase_data2 = Dict("data" => data2, "ts" => ts, "εs" => Kεs, "N" => 2, "fixedparams" => (; Δ=1), "ss_pos" => kitaev_ss)
kitaev_filename2 = savename("kitaev", kitaev_phase_data2, allowedtypes=(Number, NamedTuple), "jld2")
wsave(datadir("phase_diagram", kitaev_filename2), kitaev_phase_data2)
kitaev_phase_dataN = Dict("data" => dataN, "ts" => ts, "εs" => Kεs, "N" => N, "fixedparams" => (; Δ=1), "ss_pos" => kitaev_ss)
kitaev_filenameN = savename("kitaev", kitaev_phase_dataN, allowedtypes=(Number, NamedTuple), "jld2")
wsave(datadir("phase_diagram", kitaev_filenameN), kitaev_phase_dataN)
##
K2Data = wload(datadir("phase_diagram", "kitaev_N=2_fixedparams=(Δ = 1,).jld2"))
K40Data = wload(datadir("phase_diagram", "kitaev_N=40_fixedparams=(Δ = 1,).jld2"))
plot_data(K2Data)
plot_data(K40Data)

##
function calculate_full_phase_data(N; save, res=(100, 100), fixedparams, MaxTime=10, optimize=true, exps=range(0.1, 3, 5))
    c = FermionBdGBasis(1:N, (:↑, :↓))
    f, f!, cache = hamfunc(Hδϕ_Hε(), c, fixedparams)
    ss = find_sweet_spot((f, f!, cache), c, Hδϕ_Hε(); exps, MaxTime)
    εs = sqrt(fixedparams.Ez^2 - fixedparams.Δ^2) .+ 0.5 * range(-1, 1, length=res[1])
    δϕs = range(0, pi, length=res[2])
    iter = Iterators.product(δϕs, εs) |> collect
    mapfunc(δϕε) = fullsolve(f!(cache, δϕε), c)
    data = map(mapfunc, iter)
    join_data(ss, data, N, (εs, δϕs, ("ε", "δϕ")), "full", save)
end
function find_sweet_spot((f, f!, cache), c, optparams=Hδϕ_Hε(); exps, MaxTime, target=MPU, minexcgap=0.0)
    prob = OptProb(; hamfunc=x -> f!(cache, x), basis=c, optparams, target)
    return solve(prob, best_algs()[1]; minexcgap, maxiters=100000, MaxTime, exps)
end
function join_data(sol, data, N, (x, y, labels), prefix, save)
    phase_data = Dict("data" => data, "y" => y, "x" => x, "labels" => labels, "N" => N, "fixedparams" => fixedparams, "ss" => sol)
    if save
        filename = savename(prefix, phase_data, "jld2", allowedtypes=(Number, NamedTuple), ignores=["ss"])
        wsave(datadir("phase_diagram", "lengths", filename), phase_data)
    end
    return phase_data
end
function calculate_kitaev_phase_data(N; save, res=(100, 100))
    c = FermionBdGBasis(1:N)
    @variables t ε
    f, f! = build_function(kitaev_ham(c, ε, 1, t), [t, ε], expression=Val{false})
    cache = f([0.1, 0.1])
    ss = (; alg=:analytic, optsol=fullsolve(f([-1, 0]), c), sol=[-1, 0])
    εs = range(-3, 3, length=res[1])
    ts = range(-2, 2, length=res[2])
    iter = Iterators.product(ts, εs) |> collect
    mapfunc(tε) = fullsolve(f!(cache, tε), c)
    data = map(mapfunc, iter)
    join_data(ss, data, N, (εs, ts, ("ε", "t")), "kitaev", save)
end
##
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0.0, Ez=3)
@time data = calculate_full_phase_data(4; save=true, res=(50, 50), fixedparams, MaxTime=1, optimize=true, exps=range(0.1, 3, 5))
@time Kdata = calculate_kitaev_phase_data(2; save=false, res=(50, 50))
##
for N in 2:60
    calculate_kitaev_phase_data(N; save=true, res=(50, 50))
    calculate_full_phase_data(N; save=true, res=(50, 50), fixedparams, MaxTime=10, optimize=true, exps=range(0.1, 3, 5))
end
##
plot_LD(data)
plot_LD(Kdata)
## Full model

##
F2Data = wload(datadir("phase_diagram", "full_N=2_fixedparams=(t = 0.5, θ = QuantumDots.DiffChainParameter{Float64}(2.746801533890032), V = 0, Δ = 1, U = 0.0, Ez = 3).jld2"))
F40Data = wload(datadir("phase_diagram", "full_N=40_fixedparams=(t = 0.5, θ = QuantumDots.DiffChainParameter{Float64}(2.746801533890032), V = 0, Δ = 1, U = 0.0, Ez = 3).jld2"))
F3Data = wload(datadir("phase_diagram", "full_N=3_fixedparams=(t = 0.5, θ = QuantumDots.DiffChainParameter{Float64}(2.746801533890032), V = 0, Δ = 1, U = 0.0, Ez = 3).jld2"))
##
plot_LD(F2Data)
plot_MPU(F2Data)
plot_MP(F2Data)
plot_gap(F2Data)
##
plot_LD(F40Data)
plot_MPU(F40Data)
plot_MP(F40Data)
plot_gap(F40Data)
##

optim = Rδϕ_Rε()
sols3 = let data = data_small3, basis = FermionBdGBasis(1:3, (:↑, :↓))
    Ezs = data["Ez"]
    fixedparams = data["fixedparams"]
    optparams = data["optparams"]
    sols = Vector{Any}(undef, length(Ezs))
    for (k, (Ez, ps)) in collect(enumerate(zip(Ezs, optparams)))

        cWN = FermionBdGBasis(1:N, (:↑, :↓))
        f, f!, cache = hamfunc(optparams, cWN, fixedparams)
        prob = OptProb(; hamfunc=x -> f!(cache, x), basis=cWN, optparams, target=LD)
        sol = solve(prob, best_algs()[1]; minexcgap=0.0, maxiters=100000, MaxTime=10)

        fp = merge(fixedparams, (; Ez))
        f, f!, cache = hamfunc(optim, basis, fp)
        H = f!(cache, ps)
        sols[k] = fullsolve(H, basis)
    end
    sols
end;