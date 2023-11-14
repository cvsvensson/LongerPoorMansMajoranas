# const N = 2
# const c = FermionBasis((1, 2), (:↑, :↓); qn=QuantumDots.parity)
# hamiltonian(c; μ1, μ2, Δ, t, tratio, Ez, ϕm, ϕp, U, V) = blockdiagonal((QuantumDots.BD1_hamiltonian(c; h=Ez, t, μ=(-μ1, -μ2), Δ=Δ .* [exp(1im * ϕ / 2), exp(-1im * ϕ / 2)], Δ1=cos(ϕ / 2), θ=parameter(2atan(tratio), :diff), ϕ=0, U, V)), c)

spatial_labels(basis) = collect(unique(first.(keys(basis))))
# μsyms(N) = ntuple(i -> Symbol(:μ, i), N)
# randparams() = rand(9)
# parameternames(N) = [:t, :Δ, :V, :θ, :ϕ, :h, :U, :Δ1, μsyms(N)...]

basis(N, bdg=false) = bdg ? QuantumDots.FermionBdGBasis(1:N, (:↑, :↓)) : FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)

dictmap(f, d) = Dict(k => f(v) for (k, v) in d)
function vectorize_kwargs(c; kwargs...)
    dictmap(arg -> vectorize_arg(arg, c), kwargs)
end
vectorize_arg(arg::Number, c) = [arg for n in 1:length(spatial_labels(c))]
function vectorize_arg(arg::Vector, c)
    N = length(spatial_labels(c))
    n = length(arg)
    n == N && return arg
    2n == N && return [arg..., reverse(arg)...]
    2n - 1 == N && return [arg..., reverse(arg[begin:end-1])...]
    throw(ArgumentError("Don't know how to vectorize argument of length $n for $N sites"))
end
vectorize_arg(arg::Tuple, c) = vectorize_arg([arg...], c)

function hamiltonian(basis::FermionBdGBasis; parameters...)
    BdGMatrix(whamiltonian(basis; parameters...); check=false)
end
function hamiltonian(basis::FermionBasis; parameters...)
    blockdiagonal(whamiltonian(basis; parameters...), basis)
end

function build_whamiltonian(basis::FermionBdGBasis; parameters...)
    syms = filter(x -> x isa Num, parameters)
    bdg = BdGMatrix(whamiltonian(basis; parameters...); check=false)
    f, f! = build_function(bdg, syms, expression=Val{false})
    f2(x...) = hermitianpart!(f(x...))
    f2!(out, x...) = hermitianpart!(f!(out, x...))
    return f2, f2!
end
function build_whamiltonian(basis::FermionBasis; parameters...)
    syms = filter(x -> x isa Num, parameters)
    bdg = blockdiagonal(whamiltonian(basis; parameters...), basis)
    f, f! = build_function(bdg, syms, expression=Val{false})
    f2(x...) = hermitianpart!(f(x...))
    f2!(out, x...) = hermitianpart!(f!(out, x...))
    return f2, f2!
end

cell_labels(n, basis) = Tuple(keys(QuantumDots.cell(n, basis)))
cell_labels(basis) = Base.Fix2(cell_labels, basis)
function reduced_similarity(basis, oddvec::AbstractVector, evenvec)
    o = oddvec * oddvec'
    e = evenvec * evenvec'
    fermions = map(label -> norm(partial_trace(o - e, (label,), basis), 1), keys(basis))
    labels = cell_labels(basis)
    cells = map(n -> norm(partial_trace(o - e, labels(n), basis), 1),
        spatial_labels(basis))
    return (; fermions, cells)
end

function reduced_similarity(qps::AbstractVector{<:QuantumDots.QuasiParticle})
    basis = qps[1].basis
    N = QuantumDots.nbr_of_fermions(basis)
    labels = QuantumDots.labels(basis)
    ρeven = QuantumDots.one_particle_density_matrix(qps[1:N])
    ρodd = QuantumDots.one_particle_density_matrix(qps[vcat(1:N-1, N + 1)])
    # matrixlabels = Base.product(labels,labels)
    cell_positions(n) = [basis.position[cl] for cl in cell_labels(n, basis)]
    cinds(n) = [cell_positions(n)..., (cell_positions(n) .+ N)...]
    cell_matrices = (; even=map(n -> ρeven[cinds(n), cinds(n)], 1:div(N, 2)), odd=map(n -> ρodd[cinds(n), cinds(n)], 1:div(N, 2)))
    fermions = QuantumDots.Dictionary(labels, [norm(ρeven[[n, n + N], [n, n + N]] - ρodd[[n, n + N], [n, n + N]], 1) for n in 1:length(labels)])
    cells = QuantumDots.Dictionary(1:div(N, 2), [norm(ρeven[cinds(n), cinds(n)] - ρodd[cinds(n), cinds(n)], 1) for n in 1:div(N, 2)])
    return (; fermions, cells, cell_matrices)
end

function half_majorana_polarizations(majcoeffs, basis)
    keys1 = spatial_labels(basis)
    N = length(keys1)
    n = 1#div(N, 2)
    keys1L = keys1[1:n]
    keys1R = keys1[end:-1:end-n+1]
    keysL = filter(k -> first(k) in keys1L, keys(basis))
    keysR = filter(k -> first(k) in keys1R, keys(basis))
    left = QuantumDots.majorana_polarization(majcoeffs..., keysL)
    right = QuantumDots.majorana_polarization(majcoeffs..., keysR)
    return (; left, right)
end
function fullsolve(H, basis::FermionBasis; reduced=true, transport=missing, oddvalindex=1)
    eig = QuantumDots.diagonalize(H)
    sectors = blocks(eig)
    fullsectors = blocks(eig; full=true)
    oddvals = sectors[1].values
    evenvals = sectors[2].values
    oddvecs = fullsectors[1].vectors
    evenvecs = fullsectors[2].vectors
    oddvec = oddvecs[:, oddvalindex]
    majcoeffs = QuantumDots.majorana_coefficients(oddvec, evenvecs[:, 1], basis)
    mps = half_majorana_polarizations(majcoeffs, basis)
    reduced = reduced ? reduced_similarity(basis, oddvec, evenvecs[:, 1]) : missing
    conductance = conductance_matrix(transport, eig; basis)
    return (; gap=oddvals[oddvalindex] - first(evenvals), gapratio=gapratio(oddvals, evenvals), reduced, mps, majcoeffs, energies=(oddvals, evenvals), conductance, excgap=excgap(oddvals, evenvals))
end

function fullsolve(H, basis::FermionBdGBasis; reduced=true, transport=missing, cutoff=1e-10)
    N = QuantumDots.nbr_of_fermions(basis)
    es, ops = diagonalize(H)
    # es, ops = QuantumDots.enforce_ph_symmetry(eigen(H))
    if !QuantumDots.check_ph_symmetry(es, ops; cutoff)
        @warn "particle-hole symmetry not valid? $es \n $ops"
    end
    qps = map(op -> QuantumDots.QuasiParticle(op, basis), eachcol(ops))
    best_majorana = qps[N]
    gs_parity = QuantumDots.ground_state_parity(es, ops) |> real |> round
    gap = gs_parity == -1 ? es[N] : es[N+1]
    excgap = abs(es[N-1])
    gapratio = sign(gap)abs(gap / excgap)
    # lefthalflabels = filter(l -> Base.first(l) <= div(N, 4), keys(basis).values)
    majcoeffs = QuantumDots.majorana_coefficients(best_majorana)
    mps = half_majorana_polarizations(majcoeffs, basis)
    reduced = reduced_similarity(qps)
    return (; gap, gapratio, reduced, mps, majcoeffs, excgap, energies=es, parity=gs_parity)
end


function gapratio(oddvals, evenvals)
    δE = first(oddvals) - first(evenvals)
    Δ = min(oddvals[2], evenvals[2]) - min(first(oddvals), first(evenvals))
    return δE / Δ
end
excgap(odd, even) = min(odd[2] - odd[1], even[2] - even[1])
excgap(sol) = excgap(sol.energies...)


LD(sol) = sum(sol.reduced.cells)
LDf(sol) = sum(sol.reduced.fermions)
MP(sol) = 1 - (abs(sol.mps.left.mp) + abs(sol.mps.right.mp)) / 2
MPU(sol) = 1 - (abs(sol.mps.left.mpu) + abs(sol.mps.right.mpu)) / 2


function reflect(_p, N; pad=[])
    p = collect(_p)
    if length(p) == N
        return [p..., pad...]
    end
    if iseven(N)
        Nhalf = div(N, 2)
        if length(p) == Nhalf
            return [p..., reverse(p)..., pad...]
        else
            throw(ArgumentError("Length of p must be N or N÷2"))
        end
    end
    if isodd(N)
        Nhalf = div(N + 1, 2)
        if length(p) == Nhalf
            return [p..., reverse(p)[2:end]..., pad...]
        else
            throw(ArgumentError("Length of p must be N or N÷2"))
        end
    end
end


#
function whamiltonian_2site((c1up, c1dn), (c2up, c2dn); t, V, θϕ1, θϕ2)
    ms = hopping_rotated_nh(t, (c1up, c1dn), (c2up, c2dn), θϕ1, θϕ2)
    if iszero(V)
        return ms
    else
        return ms + V * ((numberop(c1up) + numberop(c1dn)) * (numberop(c2up) + numberop(c2dn)))
    end
end
function whamiltonian_1site((cup, cdn); ε, Ez, Δ, U)
    (ε - Ez) * numberop(cup) + (ε + Ez) * numberop(cdn) +
    pairing_nh(Δ, cup, cdn) + U * coulomb(cup, cdn)
end
pairing_nh(Δ, cup, cdn) = 2Δ * cup'cdn'
function hopping_rotated_nh(t, (c1up, c1dn), (c2up, c2dn), angles1, angles2)
    Ω = su2_rotation(angles1)' * su2_rotation(angles2)
    c1 = @SVector [c1up, c1dn]
    c2 = @SVector [c2up, c2dn]
    2t * c1' * Ω * c2
end
function whamiltonian(c; μ, Ez, t, Δ, U, V, θ)
    M = nbr_of_fermions(c)
    @assert length(cell(1, c)) == 2 "Each unit cell should have two fermions for this hamiltonian"
    N = div(M, 2)
    gv = QuantumDots.getvalue
    h1s = (whamiltonian_1site(cell(j, c); μ=gv(μ, j, N), Ez=gv(h, j, N), Δ=gv(Δ, j, N), U=gv(U, j, N)) for j in 1:N)
    h2s = (whamiltonian_2site(cell(j, c), cell(mod1(j + 1, N), c); t=gv(t, j, N; size=2), V=gv(V, j, N; size=2), θϕ1=(gv(θ, j, N), 0), θϕ2=(gv(θ, mod1(j + 1, N), N), 0)) for j in 1:N)
    return sum(h1s) + sum(h2s)
end
