abstract type AbstractOptParams end
struct Aϕ_Rε <: AbstractOptParams end
struct RΔ_Rδϕ_Rε <: AbstractOptParams end
struct Rδϕ_Rε <: AbstractOptParams end
struct Hδϕ_Hε <: AbstractOptParams end
struct Hδϕ_Aε <: AbstractOptParams end
struct Aδϕ_Aε <: AbstractOptParams end

function decompose_rϕε(ps, N=div(2length(ps), 3))
    Nhalf = div(N + 1, 2)
    Nhalf2 = div(N, 2)
    rs = ps[1:Nhalf]
    ϕs = ps[Nhalf+1:Nhalf+Nhalf2]
    εs = ps[Nhalf+Nhalf2+1:end]
    return rs, ϕs, εs
end
function decompose_ϕε(ps, N=length(ps))
    Nhalf2 = div(N, 2)
    ϕs = ps[1:Nhalf2]
    εs = ps[Nhalf2+1:end]
    return ϕs, εs
end
function decompose_allϕ_ε(ps, N)
    ϕs = ps[1:N]
    εs = ps[N+1:end]
    return ϕs, εs
end
function decompose_allδϕ_allε(ps, N)
    δϕs = ps[1:N-1]
    εs = ps[N:end]
    return δϕs, εs
end
splatter_rϕε(ps, N) = splatter_rϕε(decompose_rϕε(ps)..., N)
function splatter_rϕε(rs, δϕs, ϵs, N)
    return vcat(reflect(rs, N) .* exp.(1im .* diffreflect(δϕs, N)), ϵs)
end

hamfunc(::RΔ_Rδϕ_Rε, c, fixedparams) = hamfunc_rϕε(c, fixedparams)
hamfunc(::Rδϕ_Rε, c, fixedparams) = hamfunc_ϕε(c, fixedparams)
hamfunc(::Aδϕ_Aε, c, fixedparams) = hamfunc_allδϕ_allε(c, fixedparams)
hamfunc(::Aϕ_Rε, c, fixedparams) = hamfunc_allϕ_ε(c, fixedparams)
hamfunc(::Hδϕ_Hε, c, fixedparams) = hamfunc_Hδϕ_Hε(c, fixedparams)
hamfunc(::Hδϕ_Aε, c, fixedparams) = hamfunc_Hδϕ_Aε(c, fixedparams)
function hamfunc_rϕε(c, fixedparams)
    N = div(QuantumDots.nbr_of_fermions(c), 2)
    Nhalf = div(N + 1, 2)
    @variables Δ[1:N], εs[1:Nhalf]::Real
    params = merge(fixedparams, (; Δ=collect(Δ[1:N]), ε=reflect(εs, N)))
    f, f! = LongerPoorMansMajoranas.build_whamiltonian(c; params...)
    f2(ps) = f(splatter_rϕε(ps, N))
    f2!(out, ps) = f!(out, splatter_rϕε(ps, N))
    ps = rand(2Nhalf + div(N, 2))
    cache = get_cache(c, f2(ps))
    return f2, f2!, cache
end
function hamfunc_ϕε(c, fixedparams)
    N = div(QuantumDots.nbr_of_fermions(c), 2)
    Nhalf = div(N + 1, 2)
    Nhalf2 = div(N, 2)
    @variables Δ[1:N], εs[1:Nhalf]::Real
    params = merge(fixedparams, (; Δ=collect(Δ[1:N]), ε=reflect(εs, N)))
    f, f! = LongerPoorMansMajoranas.build_whamiltonian(c; params...)
    splatter_ϕε(ps, N) = splatter_ϕε(decompose_ϕε(ps, N)..., N)
    splatter_ϕε(δϕs, ϵs, N) = vcat(fixedparams.Δ .* exp.(1im .* diffreflect(δϕs, N)), ϵs)
    f2(ps) = f(splatter_ϕε(ps, N))
    f2!(out, ps) = f!(out, splatter_ϕε(ps, N))
    ps = rand(N)
    cache = get_cache(c, Matrix(f2(ps)))
    return f2, f2!, cache
end
function hamfunc_Hδϕ_Hε(c, fixedparams)
    N = div(QuantumDots.nbr_of_fermions(c), 2)
    @variables Δ[1:N], ε::Real
    params = merge(fixedparams, (; Δ=collect(Δ[1:N]), ε))
    f, f! = LongerPoorMansMajoranas.build_whamiltonian(c; params...)
    function splatter(δϕϵs, N)
        ϕs = (0:N-1) .* (δϕϵs[1])
        vcat(fixedparams.Δ .* exp.(1im .* ϕs), δϕϵs[2])
    end
    f2(ps) = f(splatter(ps, N))
    f2!(out, ps) = f!(out, splatter(ps, N))
    ps = rand(2)
    cache = get_cache(c, Matrix(f2(ps)))
    return f2, f2!, cache
end
function hamfunc_Hδϕ_Aε(c, fixedparams)
    N = div(QuantumDots.nbr_of_fermions(c), 2)
    @variables Δ[1:N], ε[1:N]::Real
    params = merge(fixedparams, (; Δ=collect(Δ[1:N]), ε=collect(ε[1:N])))
    f, f! = LongerPoorMansMajoranas.build_whamiltonian(c; params...)
    function splatter(δϕϵs, N)
        ϕs = (0:N-1) .* (δϕϵs[1])
        vcat(fixedparams.Δ .* exp.(1im .* ϕs), δϕϵs[2:end])
    end
    f2(ps) = f(splatter(ps, N))
    f2!(out, ps) = f!(out, splatter(ps, N))
    ps = rand(1 + N)
    cache = get_cache(c, Matrix(f2(ps)))
    return f2, f2!, cache
end
function hamfunc_allϕ_ε(c, fixedparams)
    N = div(QuantumDots.nbr_of_fermions(c), 2)
    Nhalf = div(N + 1, 2)
    # Nhalf2 = div(N, 2)
    @variables Δ[1:N], εs[1:Nhalf]::Real
    params = merge(fixedparams, (; Δ=collect(Δ[1:N]), ε=reflect(εs, N)))
    f, f! = LongerPoorMansMajoranas.build_whamiltonian(c; params...)
    splatter_allϕ_ε(ps, N) = splatter_allϕ_ε(decompose_allϕ_ε(ps, N)..., N)
    splatter_allϕ_ε(ϕs, ϵs, N) = vcat(fixedparams.Δ .* exp.(1im .* ϕs), ϵs)
    f2(ps) = f(splatter_allϕ_ε(ps, N))
    f2!(out, ps) = f!(out, splatter_allϕ_ε(ps, N))
    ps = rand(N + Nhalf)
    cache = get_cache(c, f2(ps))
    return f2, f2!, cache
end


get_initials(prob::OptProb) = get_initials(prob.optparams, div(length(prob.basis), 2))
get_ranges(prob::OptProb) = get_ranges(prob.optparams, div(length(prob.basis), 2))
initial_Aϕ(N) = zeros(N)
initial_Aδϕ(N) = zeros(N - 1)
initial_Rε(N) = zeros(div(N + 1, 2))
initial_Hε(N) = [0.0]
initial_Aε(N) = zeros(N)
initial_Hδϕ(N) = [0.0]
initial_Rδϕ(N) = 0.0 .* ones(div(N, 2))

initial_RΔ(N) = ones(div(N + 1, 2))
ranges_Aϕ(N) = [(0.0, 2.0pi) for i in 1:N]
ranges_Aδϕ(N) = [(0.0, 2.0pi) for i in 1:N-1]
ranges_Rε(N) = [100 .* (0, 1) for i in 1:div(N + 1, 2)]
ranges_Aε(N) = [100 .* (0, 1) for i in 1:N]
ranges_Rδϕ(N) = [(0.0, 1.0pi) for i in 1:div(N, 2)]
ranges_RΔ(N) = [(0.01, 10.0) for i in 1:div(N + 1, 2)]
ranges_Hε(N) = [100 .* (0, 1)]
ranges_Hδϕ(N) = [(0.0, 1.0pi)]

function get_initials(::Aϕ_Rε, N)
    vcat(initial_Aϕ(N), initial_Rε(N))
end
function get_initials(::RΔ_Rδϕ_Rε, N)
    vcat(initial_RΔ(N), initial_Rδϕ(N), initial_Rε(N))
end
function get_initials(::Hδϕ_Hε, N)
    vcat(initial_Hδϕ(N), initial_Hε(N))
end
function get_initials(::Hδϕ_Aε, N)
    vcat(initial_Hδϕ(N), initial_Aε(N))
end
function get_initials(::Rδϕ_Rε, N)
    vcat(initial_Rδϕ(N), initial_Rε(N))
end
function get_initials(::Aδϕ_Aε, N)
    vcat(initial_Aδϕ(N), initial_Aε(N))
end
function get_ranges(::Rδϕ_Rε, N)
    δϕranges = ranges_Rδϕ(N)
    εranges = ranges_Rε(N)
    vcat(δϕranges, εranges)
end
function get_ranges(::Hδϕ_Hε, N)
    δϕranges = ranges_Hδϕ(N)
    εranges = ranges_Hε(N)
    vcat(δϕranges, εranges)
end
function get_ranges(::Hδϕ_Aε, N)
    δϕranges = ranges_Hδϕ(N)
    εranges = ranges_Aε(N)
    vcat(δϕranges, εranges)
end
function get_ranges(::RΔ_Rδϕ_Rε, N)
    Δranges = ranges_RΔ(N)
    δϕranges = ranges_Rδϕ(N)
    εranges = ranges_Rε(N)
    vcat(Δranges, δϕranges, εranges)
end
function get_ranges(::Aϕ_Rε, N)
    ϕranges = ranges_Aϕ(N)
    εranges = ranges_Rε(N)
    vcat(ϕranges, εranges)
end
function get_ranges(::Aδϕ_Aε, N)
    δϕranges = ranges_Aδϕ(N)
    εranges = ranges_Aε(N)
    vcat(δϕranges, εranges)
end