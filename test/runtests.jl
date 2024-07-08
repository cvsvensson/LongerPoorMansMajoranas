using DrWatson, Test
@quickactivate :LongerPoorMansMajoranas

# Here you include files using `srcdir`
# include(srcdir("file.jl"))

# Run test suite
println("Starting tests")
ti = time()

@testset "LongerPoorMansMajoranas tests" begin end

@testset "Hamiltonians" begin
    N = 2
    δϕ = 0.3
    ϕs = (0:N-1) .* δϕ
    Δ0 = 0.4
    params = (; ε=0.1, Ez=0.2, t=0.3, Δ=Δ0 .* exp.(1im * ϕs), U=0.0, V=0.0, θ=parameter(0.7, :diff))
    fixedparams = merge(params[[:t, :U, :V, :θ, :Ez]], (; Δ=Δ0))
    c = FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
    ham1 = LongerPoorMansMajoranas.whamiltonian(c; params..., conjugate=true)
    ham2 = LongerPoorMansMajoranas.whamiltonian(c; params..., conjugate=false)
    f, f!, cache = hamfunc(Hδϕ_Hε(), c, fixedparams)
    ham3 = f([δϕ, params.ε])
    @test ham1 ≈ hermitianpart(ham2) ≈ ham3

    mbvals = eigvals(Matrix(ham1))
    mbvals = (mbvals.-mbvals[1])[2:end]

    c2 = FermionBdGBasis(1:N, (:↑, :↓))
    ham1_bdg = LongerPoorMansMajoranas.whamiltonian(c2; params..., conjugate=true)
    ham2_bdg = LongerPoorMansMajoranas.whamiltonian(c2; params..., conjugate=false)
    f_bdg, f!_bdg, cache = hamfunc(Hδϕ_Hε(), c2, fixedparams)
    ham3_bdg = f_bdg([δϕ, params.ε])
    @test ham1_bdg ≈ hermitianpart(ham2_bdg) ≈ ham3_bdg

    bdgvals = eigvals(Matrix(ham1_bdg))
    @test all(minimum(abs.(E .- mbvals)) < 1e-12 for E in bdgvals[5:8])
end

@testset "LD" begin
    N = 3
    c = FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
    cbdg = FermionBdGBasis(1:N, (:↑, :↓))
    ham = LongerPoorMansMajoranas.whamiltonian(c; ε=1, t=0.5, θ=0.7, V=0, Δ=1, U=0.0, Ez=3, conjugate=true)
    ham_bdg = LongerPoorMansMajoranas.whamiltonian(cbdg; ε=1, t=0.5, θ=0.7, V=0, Δ=1, U=0.0, Ez=3, conjugate=true)
    sol = fullsolve(ham, c)
    sol_bdg = fullsolve(ham_bdg, cbdg)
    @test sol.gap ≈ sol_bdg.gap
    @test sol.excgap ≈ sol_bdg.excgap
    @test norm(sol.majcoeffs) ≈ norm(sol_bdg.majcoeffs)
    @test LDbdg(sol) ≈ LDbdg(sol_bdg)
    @test LDf(sol) ≈ LDf(sol_bdg)
    @test sol.reduced.cells ≈ sol_bdg.reduced.cells
    @test sol.reduced.two_cells ≈ sol_bdg.reduced.two_cells
    @test collect(sol.reduced.cells_bdg) ≈ collect(sol_bdg.reduced.cells_bdg)
    @test collect(sol.reduced.two_cells_bdg) ≈ collect(sol_bdg.reduced.two_cells_bdg)
end

@testset "Equivalence" begin
    N = 3
    c = FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
    cbdg = FermionBdGBasis(1:N, (:↑, :↓))
    using Accessors
    params = (; t=0.5, θ=parameter(2atan(5), :diff), V=0.0, Δ=[1, 1, 1], U=0.0, Ez=3)
    params1 = @insert params.ε = [1, 1, 1]
    params2 = @set params.θ = parameter(-2atan(5), :diff) #This changes the sign of t/tso
    params2 = @insert params2.ε = [1, -1, 1]
    params2 = @set params2.Δ = [1, -1, 1]

    params3 = @insert params.ε = [1, -1, 1]
    params3 = @set params3.θ = parameter(-2atan(1 / 5), :diff)
    # params3 = @set params3.t = -0.5
    # params3 = @set params3.Δ = [1, -1, 1]

    # params3 = @sert params2.θ = parameter(0.7, :diff)
    ham1 = LongerPoorMansMajoranas.whamiltonian(c; params1..., conjugate=true)
    ham2 = LongerPoorMansMajoranas.whamiltonian(c; params2..., conjugate=true)
    ham3 = LongerPoorMansMajoranas.whamiltonian(c; params3..., conjugate=true)

    sol1 = fullsolve(ham1, c)
    sol2 = fullsolve(ham2, c)
    sol3 = fullsolve(ham3, c)
    # sol1.energies[1] - sol2.energies[1]
    @test LD_cells(sol1) ≈ LD_cells(sol2)
    @test !(LD_cells(sol1) ≈ LD_cells(sol3))

end

using ForwardDiff
function low_energy_projection(H)
    vals, vecs = eigen(Matrix(H))
    S = vecs[:, 1:4]
    Hermitian(S' * H * S)
end
remove_trace(x) = x - tr(x) * I / size(x, 1)
@testset "Perturbations" begin
    N = 2
    δϕ = 0.6 * pi / 2
    ϕs = (0:N-1) .* δϕ
    Δ0 = 1
    Ez = 40
    ε = 0.1 + sqrt(Ez^2 - Δ0^2)
    params = (; ε, Ez, Δ=Δ0 .* exp.(1im * ϕs), U=0.0, V=0.0, θ=parameter(2atan(5), :diff))
    c = FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
    get_ham(t; params=params) = Matrix(LongerPoorMansMajoranas.whamiltonian(c; t, params..., conjugate=true))
    a = FermionBasis(1:N; qn=QuantumDots.parity)
    get_aham_terms(t) = LongerPoorMansMajoranas.perturbative_hamiltonian_terms(a; t, δϕ=δϕ .* ones(N - 1), Δ=Δ0, params[[:ε, :Ez, :θ]]...)
    get_aham_terms(0)

    Vpert = ForwardDiff.derivative(get_ham, 1)
    Vpert2 = get_ham(1; params=(; ε=0, Ez=0, Δ=0, U=0.0, V=0.0, θ=params.θ))
    @test Vpert ≈ Vpert2

    ts = exp.(range(-10, 1, 10))
    H0 = remove_trace(get_ham(0))
    eig0 = eigen(H0)
    S = eig0.vectors
    E0 = diff(eig0.values)[1]
    # gs = S[:, 1]
    Es = [diff(eigvals(get_ham(t)))[1] for t in ts]
    # E1 = gs' * Vpert * gs

    # plot(real(diff(log.(E0 .+ E1 .* ts - Es)) ./ diff(log.(ts))))
    #  plot(real(diff(log.(abs.(E0 .- Es) ./ ts)) ./ diff(log.(ts))))

    proj = cat(I(2^N), 0I(4^N - 2^N); dims=(1, 2))
    P = S[:, 1:2^N]#S * proj * S'
    Q = S[:, 2^N+1:end]#S * (I - proj) * S'
    Heff0 = remove_trace(Hermitian(P' * H0 * P))
    aH0 = remove_trace(Hermitian(Matrix(get_aham_terms(0)[1])))
    @test eigvals(Heff0) ≈ eigvals(Matrix(aH0))

    Vd = Hermitian(P' * Vpert * P)
    Vod = P' * Vpert * Q

    aV = Hermitian(Matrix(ForwardDiff.derivative(t -> get_aham_terms(t)[2], 1)))
    @test aV ≈ get_aham_terms(1)[2]

    Ediffs0 = [norm(eigvals(Heff0) - eigvals(get_ham(t))[1:2^N] |> diff) for t in ts]
    Ediffs1 = [norm(eigvals(Heff0 + t * Vd) - eigvals(get_ham(t))[1:2^N] |> diff) for t in ts]
    plot(ts, Ediffs0, scale=:log10, markers=true)
    plot!(ts, Ediffs1, scale=:log10, markers=true)

    eiga0 = eigen(aH0)
    Sa = eiga0.vectors
    @test Sa * Heff0 * Sa' ≈ aH0

    Sa * Vd * Sa' |> norm
    aV |> norm
    # scalefactor = norm(aV) / norm(Vd)
    pretty_print(aV, a)
    pretty_print(Sa * Vd * Sa', a)
    aEdiffs0 = [norm(eigvals(Matrix(aH0)) - eigvals(get_ham(t))[1:2^N] |> diff) for t in ts]
    aEdiffs1 = [norm(eigvals(Matrix(aH0 + t * aV)) - eigvals(get_ham(t))[1:2^N] |> diff) for t in ts]
    # aEdiffs0 = [norm(eigvals(aH0) - eigvals(Matrix(sum(get_aham_terms(t))))[1:2^N] |> diff) for t in ts]
    # aEdiffs1 = [norm(eigvals(aH0 + t * aV) - eigvals(Matrix(sum(get_aham_terms(t))))[1:2^N] |> diff) for t in ts]
    plot(ts, aEdiffs0, scale=:log10, markers=true)
    plot!(ts, aEdiffs1, scale=:log10, markers=true)



    plot(real(diff(log.(abs.(Ediffs) ./ ts)) ./ diff(log.(ts))), markers=true)
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60, digits=3), " minutes")
