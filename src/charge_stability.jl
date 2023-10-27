function charge_stability_scan(parameters, dx=1, dy=1, res=100; transport=missing)
    ϵ0 = vectorize_arg(parameters.ϵ, c)
    ϵ1s = range(ϵ0[1] .- dx / 2, ϵ0[1] .+ dx / 2; length=res)
    ϵ2s = range(ϵ0[2] .- dy / 2, ϵ0[2] .+ dy / 2; length=res)
    iter = collect(Base.product(ϵ1s, ϵ2s))
    data = Folds.map((xy) -> fullsolve(hamiltonian(c; parameters..., ϵ=xy); transport), iter)
    return Dict(:data => data, :ϵs => iter, :ϵ1 => ϵ1s, :ϵ2 => ϵ2s, :parameters => parameters)
end