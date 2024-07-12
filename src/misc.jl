spatial_labels(basis) = unique(map(first, collect(keys(basis))))

get_basis(N, bdg=false) = bdg ? QuantumDots.FermionBdGBasis(1:N, (:↑, :↓)) : FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)


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

    cells_bdg = QuantumDots.Dictionary(space_labels, [norm(one_particle_density_matrix(δρ, basis, labels(n))) for n in space_labels])

    function two_cell_bdg(n)
        return norm(one_particle_density_matrix(δρ, basis, [labels(n)..., labels(n + 1)...]))
    end
    two_cells_bdg = QuantumDots.Dictionary(spatial_labels(basis)[1:end-1], [two_cell_bdg(n) for n in spatial_labels(basis)[1:end-1]])

    return (; fermions, cells, two_cells, cells_bdg, two_cells_bdg)
end

function reduced_similarity(qps::AbstractVector{<:QuantumDots.QuasiParticle}; only_cells=false)
    basis = qps[1].basis
    space_labels = spatial_labels(basis)
    N = QuantumDots.nbr_of_fermions(basis)
    labels = QuantumDots.labels(basis)
    ρeven = QuantumDots.one_particle_density_matrix(qps[1:N])
    ρodd = QuantumDots.one_particle_density_matrix(qps[vcat(1:N-1, N + 1)])
    cell_positions(n) = [basis.position[cl] for cl in cell_labels(n, basis)]
    cinds(n) = vcat(cell_positions(n), (cell_positions(n) .+ N))
    c_dot = FermionBasis(1:length(cell_labels(first(space_labels), basis)), qn=QuantumDots.parity)
    function LD_cell(n)
        inds = cinds(n)
        norm(many_body_density_matrix(ρeven[inds, inds], c_dot) - many_body_density_matrix(ρodd[inds, inds], c_dot))
    end
    cells = [LD_cell(n) for n in space_labels]
    if only_cells
        return (; cells)
    end
    fermions = QuantumDots.Dictionary(labels, [norm(ρeven[[n, n + N], [n, n + N]] - ρodd[[n, n + N], [n, n + N]]) for n in 1:length(labels)])
    function one_cell_bdg(n)
        inds = cinds(n)
        return norm(ρeven[inds, inds] - ρodd[inds, inds])
    end
    cells_bdg = QuantumDots.Dictionary(1:length(space_labels), [one_cell_bdg(n) for n in space_labels])
    function two_cell_bdg(n)
        inds = vcat(cinds(n), cinds(n + 1))
        return norm(ρeven[inds, inds] - ρodd[inds, inds])
    end
    c_2dots = FermionBasis(1:2length(cell_labels(first(space_labels), basis)))
    function LD_two_cells(n)
        inds = sort(vcat(cinds(n), cinds(n + 1)))
        norm(many_body_density_matrix(ρeven[inds, inds], c_2dots) - many_body_density_matrix(ρodd[inds, inds], c_2dots))
    end
    two_cells = [LD_two_cells(n) for n in space_labels[1:end-1]]
    two_cells_bdg = QuantumDots.Dictionary(space_labels[1:end-1], [two_cell_bdg(n) for n in space_labels[1:end-1]])
    return (; fermions, cells, two_cells, cells_bdg, two_cells_bdg)#, cell_matrices)
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

struct MBEigen{DH,B}
    eig::DH
    basis::B
end

struct BdGEigen{U,S,Q}
    values::U
    vectors::S
    qps::Q
end
function QuantumDots.diagonalize(_H::AbstractMatrix, basis::FermionBasis)
    H = Hermitian(blockdiagonal(_H, basis))
    MBEigen(diagonalize(H), basis)
end
function fullsolve(_H, basis::FermionBasis)
    eig = diagonalize(_H, basis)
    majinfo = majorana_info(eig)
    reduced = reduced_info(eig)

    return (; reduced, majinfo..., energy_info(eig)...)
end

function LD_cells(e::MBEigen)
    basis = e.basis
    fullsectors = blocks(e.eig; full=true)
    oddvec = first(eachcol(fullsectors[1].vectors))
    evenvec = first(eachcol(fullsectors[2].vectors))
    o = oddvec * oddvec'
    e = evenvec * evenvec'
    δρ = o - e
    labels = cell_labels(basis)
    space_labels = spatial_labels(basis)
    norm(norm(partial_trace(δρ, labels(n), basis)) for n in space_labels)
end
function reduced_info(eig::MBEigen)
    basis = eig.basis
    fullsectors = blocks(eig.eig; full=true)
    oddvec = first(eachcol(fullsectors[1].vectors))
    evenvec = first(eachcol(fullsectors[2].vectors))
    reduced_similarity(basis, oddvec, evenvec)
end
function majorana_info(eig::MBEigen)
    basis = eig.basis
    fullsectors = blocks(eig.eig; full=true)
    oddvec = first(eachcol(fullsectors[1].vectors))
    evenvec = first(eachcol(fullsectors[2].vectors))
    majcoeffs = QuantumDots.majorana_coefficients(oddvec, evenvec, basis)
    mps = get_majorana_polarizations(majcoeffs, basis)
    (; majcoeffs, mps)
end
function energy_info(eig::MBEigen)
    sectors = blocks(eig.eig)
    oddvals = sectors[1].values
    evenvals = sectors[2].values
    gap = first(oddvals) - first(evenvals)

    energies = (oddvals, evenvals)
    (; gap, energies, gapratio=gapratio(oddvals, evenvals), excgap=excgap(oddvals, evenvals))
end
function get_gap(eig::MBEigen)
    energy_info(eig).gap
end

function QuantumDots.diagonalize(_H::BdGMatrix, basis::FermionBdGBasis; cutoff=1e-10)
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
    return BdGEigen(es, ops, qps)
end
const two_basis = FermionBasis(1:2; qn=QuantumDots.parity)
LD_cells(e::BdGEigen) = LD_cells(e.qps)
function LD_cells(qps::AbstractVector{<:QuantumDots.QuasiParticle})
    basis = first(qps).basis
    space_labels = spatial_labels(basis)
    N = QuantumDots.nbr_of_fermions(basis)
    ρeven = QuantumDots.one_particle_density_matrix(qps[1:N])
    ρodd = QuantumDots.one_particle_density_matrix(qps[vcat(1:N-1, N + 1)])
    cell_positions(n) = [basis.position[cl] for cl in cell_labels(n, basis)]
    cinds(n) = vcat(cell_positions(n), (cell_positions(n) .+ N))
    c_dot = two_basis
    function LD_cell(n)
        inds = cinds(n)
        norm(many_body_density_matrix(ρeven[inds, inds], c_dot) - many_body_density_matrix(ρodd[inds, inds], c_dot))
    end
    norm(norm(LD_cell(n)) for n in space_labels)
end

function majorana_info(eig::BdGEigen)
    qps = eig.qps
    basis = first(qps).basis
    N = QuantumDots.nbr_of_fermions(basis)
    best_majorana = qps[N]
    majcoeffs = QuantumDots.majorana_coefficients(best_majorana)
    mps = get_majorana_polarizations(majcoeffs, basis)
    (; majcoeffs, mps)
end
function energy_info(eig::BdGEigen)
    es = eig.values
    ops = eig.vectors
    basis = first(eig.qps).basis
    N = QuantumDots.nbr_of_fermions(basis)
    gs_parity = QuantumDots.ground_state_parity(es, ops) |> real |> round
    gap = gs_parity == -1 ? es[N] : es[N+1]
    excgap = abs(es[N-1] - es[N])
    (gap, gapratio, excgap, energies=es, parity=gs_parity)
end
function get_gap(eig::BdGEigen)
    energy_info(eig).gap
end
get_gap(nt::NamedTuple) = nt.gap
function all_info(eig)
    (; majorana_info(eig)..., energy_info(eig)..., reduced=reduced_info(eig))
end
function fullsolve(_H::AbstractMatrix, basis::FermionBdGBasis)
    bdgeig = diagonalize(BdGMatrix(_H), basis)
    majinfo = majorana_info(bdgeig)
    reduced = reduced_info(bdgeig)
    return (; energy_info(bdgeig)..., majinfo..., reduced)
end
reduced_info(eig::BdGEigen) = reduced_similarity(eig.qps)
function gapratio(oddvals, evenvals)
    δE = first(oddvals) - first(evenvals)
    Δ = min(oddvals[2], evenvals[2]) - min(first(oddvals), first(evenvals))
    return δE / Δ
end
excgap(odd, even) = min(odd[2] - odd[1], even[2] - even[1])
excgap(sol) = excgap(sol.energies...)

LD_cells(sol, p=2) = norm(sol.reduced.cells, p) #maximum(sol.reduced.cells)#
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


function get_gap_gradient(eigfunc, p)
    gapfunc = get_gap ∘ eigfunc
    FiniteDiff.finite_difference_gradient(gapfunc, p)
end

function get_gap_hessian(eigfunc, p)
    gapfunc = get_gap ∘ eigfunc
    FiniteDiff.finite_difference_hessian(gapfunc, p)
end
function get_gap_derivatives(eigfunc, p)
    gapfunc = get_gap ∘ eigfunc
    gradient = FiniteDiff.finite_difference_gradient(gapfunc, p)
    hessian = FiniteDiff.finite_difference_hessian(gapfunc, p)
    (; gradient, hessian)
end