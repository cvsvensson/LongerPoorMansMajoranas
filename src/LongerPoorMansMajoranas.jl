module LongerPoorMansMajoranas
using Reexport
@reexport using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra
using SparseArrays
# using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra
using ForwardDiff
#using LinearSolve # For transport
using Folds
using StaticArrays
using Symbolics
using Optimization, OptimizationBBO, OptimizationOptimJL, OptimizationMetaheuristics
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)
using Statistics
using FiniteDiff

import AbstractDifferentiation as AD

export c, LD_cells, LDf, MP, MPU, MPI, excgap, LDbdg, MPUqd, MPqd, get_gap
export hamiltonian, cell_labels, fullsolve, reduced_similarity
export get_sweet_spot, reflect, diffreflect, OptProb, best_algs, best_alg_names
export Transport, solve
export charge_stability_scan
export Aϕ_Rε, RΔ_Rδϕ_Rε, Rδϕ_Rε, Hδϕ_Hε, Hδϕ_Aε, hamfunc
export perturbative_hamiltonian
export BestOf
export get_gap_gradient, get_gap_hessian, get_gap_derivatives, all_info
export ScheduledPenalty, ConstantPenalty, GapPenalty, MinExcGapPenalty, ScheduledOptProb

include("misc.jl")
include("optimize.jl")
include("hamfunc.jl")
include("charge_stability.jl")
include("perturbation_julia_expressions.jl")
include("hamiltonians.jl")
end