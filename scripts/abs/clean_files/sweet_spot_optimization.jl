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
synceddir(args...) = joinpath(ENV["Dropbox"], "data", "LongerPoorMans", args...)
includet(scriptsdir("abs", "phase_plots.jl"))
includet(scriptsdir("abs", "phase_misc.jl"))

##
data = []
# initials = [[1.9, 2.88], [2.2, 2.95]]
initials = [[1.9, 2.9], [1, 1]]
fixedparams = (; t=0.5, θ=parameter(2atan(5), :diff), V=0, Δ=1, U=0.0, Ez=3)
for N in 2:20
    alg = BestOf(best_algs()[1:end-1])
    target = LDbdg
    #calculate_kitaev_phase_data(N; save=true, res=(50, 50))
    d = calculate_full_phase_data(N; bdg=true, save=true, res=(2, 2), fixedparams, MaxTime=10sqrt(N), target, optimize=true, exps=[-10], folder="ss_bdg_noexcgap_bestof_nodeg", final_NM=false, minexcgap=0, alg, initials=initials[1])
    # initials[1] = initials[2]
    # initials[2] = d["ss"].sol
    push!(data, d)
end
##
for N in [20, 40]
    #calculate_kitaev_phase_data(N; save=true, res=(50, 50))
    d = calculate_full_phase_data(N; bdg=true, save=true, res=(500, 500), fixedparams, optimize=false, folder="high_res")
end
##
fig = plot(; yscale=:log);
let
    f = d -> LDf(d.optsol)
    f = d -> LDbdgmax(d.optsol)
    ss = map(d -> d["ss"], data)
    f = d -> abs(d.optsol.gap)
    # f = d -> d.sol[1]
    ds = map(x -> f.(x["ss"].all_ss), data)
    d = map(f, ss)
    println(sum(d))
    plot!(fig, d, ls=:solid, ylims=(0.0000001, 1))
    dsp = [[(i, d) for d in ds[i]] for i in eachindex(ds)]
    println(dsp)
    scatter!(fig, vcat(dsp...))

end
fig
