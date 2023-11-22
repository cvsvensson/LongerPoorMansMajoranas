module LongerPoorMansMajoranas
using Reexport
@reexport using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra
# using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra
using ForwardDiff, LinearSolve # For transport
using Folds
using StaticArrays
using Symbolics
using Optimization, OptimizationBBO, OptimizationOptimJL, OptimizationMetaheuristics

export c, LD, LDf, MP, MPU, excgap
export hamiltonian, cell_labels, fullsolve, reduced_similarity
export get_sweet_spot, reflect, diffreflect, OptProb, best_algs, best_alg_names
export Transport, solve
export charge_stability_scan

include("misc.jl")
include("optimize.jl")
include("transport.jl")
include("charge_stability.jl")
end