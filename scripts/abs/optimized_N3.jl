sols = let fixedparams, c, εs, δϕs, iter, mapfunc, data, ss, f, fK
    δϕ = [0.754094, -0.0226959 - 0.754094]
    ϕ = push!(cumsum(δϕ), 0)
    Δ0 = 1
    t = 0.0879
    Ez = 1.3143
    θ = parameter(2atan(0.995962), :diff)
    fixedparams = (; t, θ, V=0, Δ=Δ0 .* exp.(1im * ϕ), U=0, Ez)

    c = FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
    a = FermionBasis(1:N; qn=QuantumDots.parity)
    f(ε) = blockdiagonal(LongerPoorMansMajoranas.whamiltonian(c; ε, fixedparams..., conjugate=true), c)
    # fK(ε) = perturbative_hamiltonian_terms(c; ε, fixedparams...)

    εs = sqrt(fixedparams.Ez^2 - Δ0^2) .+ 0.7 * fixedparams.t * range(-1, 1, length=100)

    fixedparams2 = (; δϕ, t, θ, V=0, Δ=Δ0, U=0, Ez)
    psolsall = []
    for ε in εs
        fp = filter(kv -> kv[1] ∉ (:U, :V), pairs(fixedparams2)) |> NamedTuple
        params = merge(fp, (; ε))
        HKs = LongerPoorMansMajoranas.perturbative_hamiltonian_terms(a; params...) .|> Hermitian
        # display(diagonalize(sum(HKs[1:2])))
        psols = [fullsolve(blockdiagonal(sum(HKs[1:k]), a), a) for k in 1:3]
        push!(psolsall, psols)
    end
    println(fixedparams)

    mapfunc(δϕε) = fullsolve(f(δϕε), c)
    data = map(mapfunc, εs)
    hcat(stack(psolsall) |> permutedims, data)
end
##
plot(map(LD, sols))
plot(map(x -> x.gap, sols))
plot(map(MPU, sols))
