using DrWatson, Test
@quickactivate :LongerPoorMansMajoranas
using Accessors

@testset "Consistency between bdg and many-body" begin
    N = 3
    c = FermionBasis(1:N, (:↑, :↓); qn=QuantumDots.parity)
    cbdg = FermionBdGBasis(1:N, (:↑, :↓))
    fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0.0, Ez=3)
    f, f!, cache = hamfunc(Hδϕ_Hε(), c, fixedparams)
    fbdg, fbdg!, cachebdg = hamfunc(Hδϕ_Hε(), cbdg, fixedparams)

    @test f([0.1, 0.2]) ≈ f!(cache, [0.1, 0.2])
    @test fbdg([0.1, 0.2]) ≈ fbdg!(cachebdg, [0.1, 0.2])

    sol = fullsolve(cache, c)
    solbdg = fullsolve(cachebdg, cbdg)

    @test collect(sol.reduced.fermions) ≈ collect(solbdg.reduced.fermions)
    @test collect(sol.reduced.cells_bdg) ≈ collect(solbdg.reduced.cells_bdg)
    @test (sol.mps.dots[1].mpu ≈ solbdg.mps.dots[1].mpu)
    @test collect(sol.reduced.two_cells_bdg) ≈ collect(solbdg.reduced.two_cells_bdg)
end

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
