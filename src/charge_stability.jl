function charge_stability_scan(parameters, dx=1, dy=1, res=100; transport=missing, basis)
    μ0 = vectorize_arg(parameters.μ, basis)
    μ1s = range(μ0[1] .- dx / 2, μ0[1] .+ dx / 2; length=res)
    μ2s = range(μ0[2] .- dy / 2, μ0[2] .+ dy / 2; length=res)
    iter = collect(Base.product(μ1s, μ2s))
    data = Folds.map((xy) -> fullsolve(hamiltonian(basis; parameters..., μ=xy), basis; transport), iter)
    return Dict(:data => data, :μs => iter, :μ1 => μ1s, :μ2 => μ2s, :parameters => parameters)
end