
## Transport
struct Transport{T,NT}
    type::T
    parameters::NT
end


QuantumDots.conductance_matrix(::Missing, eig, basis; kwargs...) = missing
function QuantumDots.conductance_matrix(t::Transport, eig, basis)
    leads = get_leads(basis, t.parameters...)
    # system = #t.type(QuantumDots.diagonalize(QuantumDots.OpenSystem(eig, leads)))
    system = t.type(eig, leads)
    QuantumDots.conductance_matrix(AD.FiniteDifferencesBackend(FiniteDifferences.central_fdm(2, 1)), system)
end
# function get_leads(c, T, μ, Γ=1)
#     l = spatial_labels(c)
#     left = QuantumDots.CombinedLead((Γ * c[l[1], :↑]', Γ * c[l[1], :↓]'); T, μ=μ[1])
#     right = QuantumDots.CombinedLead((Γ * c[l[end], :↑]', Γ * c[l[end], :↓]'); T, μ=μ[2])
#     return (; left, right)
# end
function get_leads(c, T, μ, Γ=1)
    l = spatial_labels(c)
    leftfs = Γ .* Tuple(QuantumDots.cell(l[1], c))
    rightfs = Γ .* Tuple(QuantumDots.cell(l[end], c))
    left = QuantumDots.CombinedLead(leftfs; T, μ=μ[1])
    right = QuantumDots.CombinedLead(rightfs; T, μ=μ[2])
    return (; left, right)
end

