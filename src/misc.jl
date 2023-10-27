# const N = 2
# const c = FermionBasis((1, 2), (:↑, :↓); qn=QuantumDots.parity)
# hamiltonian(c; ϵ1, ϵ2, Δ, t, tratio, Ez, ϕm, ϕp, U, V) = blockdiagonal((QuantumDots.BD1_hamiltonian(c; h=Ez, t, μ=(-ϵ1, -ϵ2), Δ=Δ .* [exp(1im * ϕ / 2), exp(-1im * ϕ / 2)], Δ1=cos(ϕ / 2), θ=parameter(2atan(tratio), :diff), ϕ=0, U, V)), c)

spatial_labels(basis) = collect(unique(first.(keys(basis))))
μsyms(N) = ntuple(i -> Symbol(:μ, i), N)
randparams() = rand(9)
parameternames(N) = [:t, :Δ, :V, :θ, :ϕ, :h, :U, :Δ1, μsyms(N)...]

basis(N, bdg=false) = bdg ? QuantumDots.FermionBdGBasis(1:N, (:↑, :↓)) : FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)

dictmap(f, d) = Dict(k => f(v) for (k, v) in d)
function vectorize_kwargs(c; kwargs...)
    dictmap(arg -> vectorize_arg(arg, c), kwargs)
end
vectorize_arg(arg::Number, c) = [arg for n in 1:length(spatial_labels(c))]
vectorize_arg(arg::Vector, c) = arg
vectorize_arg(arg::Tuple, c) = arg

function hamiltonian(c; kwargs...)
    hamiltonian_vecargs(c; vectorize_kwargs(c; kwargs...)...)
end
const η = (; ↑=1, ↓=-1)
function hamiltonian_vecargs(c; ϵ, Δ, t, Ez, ϕm, ϕp, U, θ)
    sico = sincos.(θ / 2)
    si = first.(sico)
    co = last.(sico)
    ep = exp.(1im * ϕp / 2)
    cϕ = cos.(ϕm / 2)
    bar(σ) = σ == :↑ ? :↓ : :↑
    N = length(spatial_labels(c))
    H = sum(ϵ[k] * c[k, σ]'c[k, σ] for (k, σ) in keys(c)) +
        sum(U[k] * c[k, :↑]'c[k, :↑] * c[k, :↓]'c[k, :↓] for k in spatial_labels(c)) +
        sum(-Ez[k] * η[σ] * c[k, σ]'c[k, σ] for (k, σ) in keys(c)) +
        sum(t[k] * co[k] * c[k, σ]'c[k+1, σ] + η[σ] * t[k] * si[k] * c[k, σ]'c[k+1, bar(σ)] + hc for (k, σ) in keys(c) if k < N) +
        sum(Δ[k] * ep[k] * cϕ[k] * c[k, :↑]'c[k, :↓]' + hc for k in 1:N) +
        sum(η[σ]Δ[k] * ep[k] * cϕ[k] * (η[bar(σ)] * si[k] * c[k, σ]'c[k+1, σ]' + co[k] * c[k, σ]'c[k+1, bar(σ)]') + hc for (k, σ) in keys(c) if k < N)
    blockdiagonal(Hermitian(H), c)
end


cell_labels(n, basis) = Tuple(keys(QuantumDots.cell(n, basis)))
cell_labels(basis) = Base.Fix2(cell_labels, basis)
function reduced_similarity(basis, oddvec::AbstractVector{T}, evenvec) where {T}
    o = oddvec * oddvec'
    e = evenvec * evenvec'
    fermions = map(label -> norm(partial_trace(o - e, (label,), basis), 1), keys(basis))
    labels = cell_labels(basis)
    cells = map(n -> norm(partial_trace(o - e, labels(n), basis), 1),
        spatial_labels(basis))
    return (; fermions, cells)
end

function half_majorana_polarizations(majcoeffs, basis)
    keys1 = spatial_labels(basis)
    N = length(keys1)
    n = div(N, 2)
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
    return (; gap=oddvals[oddvalindex] - first(evenvals), gapratio=gapratio(oddvals, evenvals), reduced, mps, majcoeffs, energies=(oddvals, evenvals), conductance)
end

function fullsolve(H, basis::FermionBdGBasis; reduced=true, transport=missing, oddvalindex=1, cutoff = 1e-10)
    N = QuantumDots.nbr_of_fermions(basis)
    es, ops = QuantumDots.enforce_ph_symmetry(eigen(H))
    if !QuantumDots.check_ph_symmetry(es, ops; cutoff)
        @warn "particle-hole symmetry not valid? $es \n $ops"
    end
    qps = map(op -> QuantumDots.QuasiParticle(op, basis), eachcol(ops))
    best_majorana = qps[N]
    gs_parity = QuantumDots.ground_state_parity(es, ops)
    gap = gs_parity == -1 ? es[N] : es[N+1]
    gapratio = sign(gap)abs(gap / (es[N-1]))
    lefthalflabels = filter(l -> Base.first(l) <= div(N, 4), keys(basis).values)
    majcoeffs = QuantumDots.majorana_wavefunctions(best_majorana)
    md1, md2 = QuantumDots.majorana_densities(best_majorana, lefthalflabels)
    mpu = md1 - md2
    mp = mpu / (md1 + md2)
    reduced = reduced_similarity(qps)
    return (; gap, gapratio, reduced, mp, mpu, majcoeffs, energies=es)

end


function gapratio(oddvals, evenvals)
    δE = first(oddvals) - first(evenvals)
    Δ = min(oddvals[2], evenvals[2]) - min(first(oddvals), first(evenvals))
    return δE / Δ
end
excgap(odd, even) = even[2] - even[1]#min(odd[2] - odd[1], even[2] - even[1])
excgap(sol) = excgap(sol.energies...)


LD(sol) = sum(sol.reduced.cells)
LDf(sol) = sum(sol.reduced.fermions)
MP(sol) = 1 - (abs(sol.mps.left.mp) + abs(sol.mps.right.mp)) / 2
MPU(sol) = 1 - (abs(sol.mps.left.mpu) + abs(sol.mps.right.mpu)) / 2

