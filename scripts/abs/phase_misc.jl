kitaev_ham(c, ε, Δ, t) = BdGMatrix(QuantumDots.kitaev_hamiltonian(c; μ=-ε, t, Δ); check=false)
function calculate_full_phase_data(N; save, res=(100, 100), fixedparams, MaxTime=10, optimize=true, exps=range(0.1, 3, 5), folder, scale = 1)
    c = FermionBasis(1:N, (:↑, :↓); qn = QuantumDots.parity)
    # c = FermionBdGBasis(1:N, (:↑, :↓))
    f, f!, cache = hamfunc(Hδϕ_Hε(), c, fixedparams)
    ss = optimize ? find_sweet_spot((f, f!, cache), c, Hδϕ_Hε(); exps, MaxTime) : missing
    εs = sqrt(fixedparams.Ez^2 - fixedparams.Δ^2) .+ scale * 0.7 * fixedparams.t * range(-1, 1, length=res[1])
    δϕs = range(0, pi, length=res[2])
    iter = Iterators.product(δϕs, εs) |> collect
    mapfunc(δϕε) = fullsolve(f!(cache, δϕε), c)
    data = map(mapfunc, iter)
    join_data(ss, data, N, (εs, δϕs, ("ε", "δϕ")), "full", save, folder)
end
function find_sweet_spot(N; MaxTime, exps=range(0.1, 3, 5))
    c = FermionBdGBasis(1:N, (:↑, :↓))
    f, f!, cache = hamfunc(Hδϕ_Hε(), c, fixedparams)
    find_sweet_spot((f, f!, cache), c, Hδϕ_Hε(); exps, MaxTime)
end
function find_sweet_spot((f, f!, cache), c, optparams=Hδϕ_Hε(); exps, MaxTime, target=MPU, minexcgap=0.0)
    prob = OptProb(; hamfunc=x -> f!(cache, x), basis=c, optparams, target)
    return solve(prob, best_algs()[1]; minexcgap, maxiters=100000, MaxTime, exps)
end
function join_data(sol, data, N, (x, y, labels), prefix, save, folder)
    phase_data = Dict("data" => data, "y" => y, "x" => x, "labels" => labels, "N" => N, "fixedparams" => fixedparams, "ss" => sol)
    if save
        filename = savename(prefix, phase_data, "jld2", allowedtypes=(Number, NamedTuple), ignores=["ss"])
        wsave(datadir("phase_diagram", folder, filename), phase_data)
    end
    return phase_data
end
function calculate_kitaev_phase_data(N; save, res=(100, 100), folder)
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
    join_data(ss, data, N, (εs, ts, ("ε", "t")), "kitaev", save, folder)
end