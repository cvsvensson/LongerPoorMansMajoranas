spatial_labels(basis) = collect(unique(first.(keys(basis))))

basis(N, bdg=false) = bdg ? QuantumDots.FermionBdGBasis(1:N, (:↑, :↓)) : FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)


function get_symlist(parameters)
    syms = filter(x -> x[2] isa Num || x[2] isa Vector{Num}, parameters) |> values |> collect
    unique(reduce(vcat, syms))
end
function build_whamiltonian(basis::FermionBdGBasis; parameters...)
    bdg = Matrix(whamiltonian(basis; conjugate=false, parameters...))
    symlist = get_symlist(parameters)
    f, f! = build_function(bdg, symlist, expression=Val{false})
    f2(x...) = hermitianpart(f(x...)) |> BdGMatrix
    f2!(out, x...) = (f!(out, x...); hermitianpart!(out) |> BdGMatrix)
    return f2, f2!
end
function build_whamiltonian(basis::FermionBasis; parameters...)
    symlist = get_symlist(parameters)
    bd = blockdiagonal(Matrix(whamiltonian(basis; conjugate=false, parameters...)), basis)
    f, f! = build_function(bd, symlist, expression=Val{false})
    f2(x...) = hermitianpart!(f(x...))
    f2!(out, x...) = (f!(out, x...); hermitianpart!(out))
    return f2, f2!
end

cell_labels(n, basis) = Tuple(keys(QuantumDots.cell(n, basis)))
cell_labels(basis) = Base.Fix2(cell_labels, basis)
function reduced_similarity(basis, oddvec::AbstractVector, evenvec)
    o = oddvec * oddvec'
    e = evenvec * evenvec'
    δρ = o - e
    fermions = map(label -> norm(partial_trace(δρ, (label,), basis)), keys(basis))
    labels = cell_labels(basis)
    space_labels = spatial_labels(basis)

    cells = map(n -> norm(partial_trace(δρ, labels(n), basis), 2), space_labels)
    two_cells = map((n1, n2) -> norm(partial_trace(δρ, sort(union(labels(n1), labels(n2)); by=last), basis), 2), space_labels, Iterators.drop(space_labels, 1))

    local_fermions(n) = reduce(vcat, [basis[l], basis[l]'] for l in labels(n))
    # @time cells_bdg = QuantumDots.Dictionary(space_labels, [sqrt(sum(abs2, tr(δρ * (f1 * f2)) for (f1, f2) in Base.product(local_fermions(n), local_fermions(n)))) for n in space_labels])
    cells_bdg = QuantumDots.Dictionary(space_labels, [norm(one_particle_density_matrix(δρ, basis, labels(n))) for n in space_labels])

    fermions_opdm = QuantumDots.Dictionary(keys(basis), [norm([tr(δρ * (f' * f)), tr(δρ * (f * f'))]) for f in basis])

    function two_cell_bdg(n)
        return norm(one_particle_density_matrix(δρ, basis, [labels(n)..., labels(n + 1)...]))
    end
    two_cells_bdg = QuantumDots.Dictionary(spatial_labels(basis)[1:end-1], [two_cell_bdg(n) for n in spatial_labels(basis)[1:end-1]])

    return (; fermions, cells, two_cells, cells_bdg, two_cells_bdg, fermions_opdm)
end

function reduced_similarity(qps::AbstractVector{<:QuantumDots.QuasiParticle})
    basis = qps[1].basis
    space_labels = spatial_labels(basis)
    N = QuantumDots.nbr_of_fermions(basis)
    labels = QuantumDots.labels(basis)
    ρeven = QuantumDots.one_particle_density_matrix(qps[1:N])
    ρodd = QuantumDots.one_particle_density_matrix(qps[vcat(1:N-1, N + 1)])
    cell_positions(n) = [basis.position[cl] for cl in cell_labels(n, basis)]
    cinds(n) = vcat(cell_positions(n), (cell_positions(n) .+ N))
    fermions = QuantumDots.Dictionary(labels, [norm(ρeven[[n, n + N], [n, n + N]] - ρodd[[n, n + N], [n, n + N]]) for n in 1:length(labels)])

    function one_cell_bdg(n)
        inds = cinds(n)
        return norm(ρeven[inds, inds] - ρodd[inds, inds])
    end
    cells_bdg = QuantumDots.Dictionary(1:div(N, 2), [one_cell_bdg(n) for n in space_labels])
    function two_cell_bdg(n)
        inds = vcat(cinds(n), cinds(n + 1))
        return norm(ρeven[inds, inds] - ρodd[inds, inds])
    end
    c_dot = FermionBasis(1:length(cell_labels(first(space_labels), basis)), qn=QuantumDots.parity)
    function LD_cell(n)
        inds = cinds(n)
        norm(many_body_density_matrix(ρeven[inds, inds], c_dot) - many_body_density_matrix(ρodd[inds, inds], c_dot))
    end
    c_2dots = FermionBasis(1:2length(cell_labels(first(space_labels), basis)))
    function LD_two_cells(n)
        inds = sort(vcat(cinds(n), cinds(n + 1)))
        norm(many_body_density_matrix(ρeven[inds, inds], c_2dots) - many_body_density_matrix(ρodd[inds, inds], c_2dots))
    end
    cells = [LD_cell(n) for n in space_labels]
    two_cells = [LD_two_cells(n) for n in space_labels[1:end-1]]
    two_cells_bdg = QuantumDots.Dictionary(space_labels[1:end-1], [two_cell_bdg(n) for n in space_labels[1:end-1]])
    return (; fermions, cells_bdg, two_cells_bdg, cells, two_cells)#, cell_matrices)
end

function get_majorana_polarizations(majcoeffs, basis)
    keys1 = spatial_labels(basis)
    n = 1#div(N, 2)
    keys1L = keys1[1:n]
    keys1R = keys1[end:-1:end-n+1]
    keysL = filter(k -> first(k) in keys1L, keys(basis))
    keysR = filter(k -> first(k) in keys1R, keys(basis))
    left = QuantumDots.majorana_polarization(majcoeffs..., keysL)
    right = QuantumDots.majorana_polarization(majcoeffs..., keysR)
    dots = [QuantumDots.majorana_polarization(majcoeffs..., filter(k -> first(k) == key, keys(basis))) for key in keys1]
    return (; left, right, dots)
end
function fullsolve(_H, basis::FermionBasis; reduced=true, oddvalindex=1)
    H = Hermitian(blockdiagonal(_H, basis))
    eig = QuantumDots.diagonalize(H)
    sectors = blocks(eig)
    fullsectors = blocks(eig; full=true)
    oddvals = sectors[1].values
    evenvals = sectors[2].values
    oddvecs = fullsectors[1].vectors
    evenvecs = fullsectors[2].vectors
    oddvec = oddvecs[:, oddvalindex]
    majcoeffs = QuantumDots.majorana_coefficients(oddvec, evenvecs[:, 1], basis)
    mps = get_majorana_polarizations(majcoeffs, basis)
    reduced = reduced ? reduced_similarity(basis, oddvec, evenvecs[:, 1]) : missing
    return (; gap=oddvals[oddvalindex] - first(evenvals), gapratio=gapratio(oddvals, evenvals), reduced, mps, majcoeffs, energies=(oddvals, evenvals), excgap=excgap(oddvals, evenvals))
end
using SparseArrays
function fullsolve(_H::AbstractMatrix, basis::FermionBdGBasis; reduced=true, cutoff=1e-10)
    H = BdGMatrix(_H |> Hermitian; check=false)
    N = QuantumDots.nbr_of_fermions(basis)
    es, ops = try
        diagonalize(H)
    catch y
        println(sparse(QuantumDots.bdg_to_skew(H)))
        rethrow()
    end
    if !QuantumDots.check_ph_symmetry(es, ops; cutoff)
        @warn "particle-hole symmetry not valid?" #$es \n $ops"
        @debug "particle-hole symmetry not valid?" es ops inds = Iterators.take(eachindex(es), N) p = sortperm(es, by=QuantumDots.energysort)
        # @debug "$(sum(abs(es[p[i]] + es[p[QuantumDots.quasiparticle_adjoint_index(i, N)]]) for i in inds))"
        # @debug "$(norm(ops' * ops - I)))"
    end
    qps = map(op -> QuantumDots.QuasiParticle(op, basis), eachcol(ops))
    best_majorana = qps[N]
    gs_parity = QuantumDots.ground_state_parity(es, ops) |> real |> round
    gap = gs_parity == -1 ? es[N] : es[N+1]
    excgap = abs(es[N-1] - es[N])
    gapratio = sign(gap)abs(gap / excgap)
    majcoeffs = QuantumDots.majorana_coefficients(best_majorana)
    mps = get_majorana_polarizations(majcoeffs, basis)
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

LD(sol, p=2) = norm(sol.reduced.cells, p) #maximum(sol.reduced.cells)#
LDmax(sol) = maximum(sol.reduced.cells)
LDf(sol, p=2) = norm(sol.reduced.fermions, p)
LDfmax(sol) = maximum(sol.reduced.fermions)
LDbdgmax(sol) = maximum(sol.reduced.cells_bdg)
LDbdg(sol, p=2) = norm(sol.reduced.cells_bdg, p)
MP(sol) = 1 - (abs(sol.mps.left.mp) + abs(sol.mps.right.mp)) / 2
MPqd(sol) = 1 - mean(abs ∘ (x -> x.mp), sol.mps.dots)
MPUqd(sol) = 1 - mean(abs ∘ (x -> x.mpu), sol.mps.dots)
MPU(sol) = 1 - (abs(sol.mps.left.mpu) + abs(sol.mps.right.mpu)) / 2


function diffreflect(_dp, N; p0=zero(eltype(_dp)))
    p = collect(_dp)
    if N == 2
        @assert length(p) == 1
        return [p0, p...]
    end
    @assert length(p) == div(N, 2) "$(length(p)) !== $(div(N, 2))"
    p = isodd(N) ? accumulate(+, [p0, p..., reverse(p)...]) : accumulate(+, [p0, p..., reverse(p[1:end-1])...])
    return p
end
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


function get_gap_gradient(fs, p)
    gapfunc = (x -> x.gap) ∘ fs
    FiniteDiff.finite_difference_gradient(gapfunc, p)
end

function get_gap_hessian(fs, p)
    gapfunc = (x -> x.gap) ∘ fs
    FiniteDiff.finite_difference_hessian(gapfunc, p)
end
function get_gap_derivatives(fs, p)
    gapfunc = (x -> x.gap) ∘ fs
    gradient = FiniteDiff.finite_difference_gradient(gapfunc, p)
    hessian = FiniteDiff.finite_difference_hessian(gapfunc, p)
    (; gradient, hessian)
end