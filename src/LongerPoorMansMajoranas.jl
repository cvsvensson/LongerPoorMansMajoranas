module LongerPoorMansMajoranas
using Reexport
@reexport using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra
# using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra
using ForwardDiff, LinearSolve # For transport
using Folds
using StaticArrays
using Symbolics
using Optimization, OptimizationBBO, OptimizationOptimJL, OptimizationMetaheuristics
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)
using Statistics

import AbstractDifferentiation as AD
using FiniteDifferences

export c, LD, LDf, MP, MPU, MPI, excgap
export hamiltonian, cell_labels, fullsolve, reduced_similarity
export get_sweet_spot, reflect, diffreflect, OptProb, best_algs, best_alg_names
export Transport, solve
export charge_stability_scan
export Aϕ_Rε, RΔ_Rδϕ_Rε, Rδϕ_Rε, Hδϕ_Hε, Hδϕ_Aε, hamfunc
export perturbative_hamiltonian
export BestOf

include("misc.jl")
include("optimize.jl")
include("transport.jl")
include("charge_stability.jl")
include("perturbation_julia_expressions.jl")
include("hamiltonians.jl")
end