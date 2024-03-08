kitaev_ham(c, ε, Δ, t) = BdGMatrix(QuantumDots.kitaev_hamiltonian(c; μ=-ε, t, Δ); check=false, kwargs...)
function calculate_full_phase_data(N; bdg, save, res=(100, 100), fixedparams, MaxTime=10, optimize=true, exps=range(0.1, 3, 5), folder, scale=1, transport=missing, kwargs...)
    # c = FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
    c = bdg ? FermionBdGBasis(1:N, (:↑, :↓)) : FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
    f, f!, cache = hamfunc(Hδϕ_Hε(), c, fixedparams)
    ss = optimize ? find_sweet_spot((f, f!, cache), c, Hδϕ_Hε(); exps, MaxTime, kwargs...) : missing
    εs = sqrt(fixedparams.Ez^2 - fixedparams.Δ^2) .+ scale * 0.7 * fixedparams.t * range(-1, 1, length=res[1])
    δϕs = range(0, pi, length=res[2])
    iter = Iterators.product(δϕs, εs) |> collect
    mapfunc(δϕε) = fullsolve(f!(cache, δϕε), c; transport)
    data = map(mapfunc, iter)
    join_data(ss, data, N, (εs, δϕs, ("ε", "δϕ")), "full", save, folder)
end

function calculate_refl_phase_data(N, εmid; save, res=(100, 100), fixedparams, MaxTime=10, optimize=true, exps=range(0.1, 3, 5), folder, scale=1, transport=missing)
    # c = FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
    c = FermionBdGBasis(1:N, (:↑, :↓))
    f, f!, cache = hamfunc(Rδϕ_Rε(), c, fixedparams)
    ss = optimize ? find_sweet_spot((f, f!, cache), c, Rδϕ_Rε(); exps, MaxTime) : missing
    εs = sqrt(fixedparams.Ez^2 - fixedparams.Δ^2) .+ scale * 0.7 * fixedparams.t * range(-1, 1, length=res[1])
    δϕs = range(0, pi, length=res[2])
    iter = Iterators.product(δϕs, εs) |> collect
    mapfunc(δϕε) = fullsolve(f!(cache, [δϕε[1], δϕε[2], εmid]), c; transport)
    data = map(mapfunc, iter)
    join_data(ss, data, N, (εs, δϕs, ("ε", "δϕ")), "full", save, folder)
end

function conductance_sweep(fixedparams, ss, εs, Vs, T)
    δϕ = ss[1]
    ε0 = ss[2]
    c = FermionBasis(1:3, (:↑, :↓), qn=QuantumDots.parity)
    f, f!, cache = hamfunc(Hδϕ_Aε(), c, fixedparams)
    tl(V) = Transport(QuantumDots.PauliSystem, (; T, μ=(V, 0.0)))
    tr(V) = Transport(QuantumDots.PauliSystem, (; T, μ=(0.0, V)))
    iter = Iterators.product(εs, Vs)
    onel = map((εV) -> fullsolve(f!(cache, [δϕ, εV[1], ε0, ε0]), c; transport=tl(εV[2])), iter)
    twol = map((εV) -> fullsolve(f!(cache, [δϕ, εV[1], εV[1], ε0]), c; transport=tl(εV[2])), iter)
    twor = map((εV) -> fullsolve(f!(cache, [δϕ, εV[1], εV[1], ε0]), c; transport=tr(εV[2])), iter)
    threel = map((εV) -> fullsolve(f!(cache, [δϕ, εV[1], εV[1], εV[1]]), c; transport=tl(εV[2])), iter)
    return (; onel, twol, twor, threel)
end
function pert_conductance_sweep(fixedparams, ss, εs, Vs, T)
    a = FermionBasis(1:3, qn=QuantumDots.parity)
    δϕ = ones(2) .* ss[1]
    ε0 = ss[2]
    tl(V) = Transport(QuantumDots.PauliSystem, (; T, μ=(V, 0.0)))
    tr(V) = Transport(QuantumDots.PauliSystem, (; T, μ=(0.0, V)))
    iter = Iterators.product(εs, Vs)
    fp = filter(kv -> kv[1] ∉ (:U, :V), pairs(fixedparams)) |> NamedTuple

    function get_nt(M)
        onel = map((εV) -> fullsolve(perturbative_hamiltonian(a, M; fp..., δϕ, ε=[εV[1], ε0, ε0]), a; transport=tl(εV[2])), iter)
        twol = map((εV) -> fullsolve(perturbative_hamiltonian(a, M; fp..., δϕ, ε=[εV[1], εV[1], ε0]), a; transport=tl(εV[2])), iter)
        twor = map((εV) -> fullsolve(perturbative_hamiltonian(a, M; fp..., δϕ, ε=[εV[1], εV[1], ε0]), a; transport=tr(εV[2])), iter)
        threel = map((εV) -> fullsolve(perturbative_hamiltonian(a, M; fp..., δϕ, ε=[εV[1], εV[1], εV[2]]), a; transport=tl(εV[2])), iter)
        return (; onel, twol, twor, threel)
    end
    map(get_nt, 1:2)
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