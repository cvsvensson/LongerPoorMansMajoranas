using DrWatson
@quickactivate :LongerPoorMansMajoranas
using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra#, BlackBoxOptim
using Plots
using Symbolics
using Folds
using Accessors
using JLD2
using DataFrames
using LaTeXStrings
using Test
## Consistency between bdg and many-body
@testset "Consistency between bdg and many-body" begin
    N = 2
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
    @test (sol.mps.singles[1].mpu ≈ solbdg.mps.singles[1].mpu)
end