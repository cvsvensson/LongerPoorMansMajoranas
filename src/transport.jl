
## Transport
struct Transport{T,NT}
    type::T
    parameters::NT
end
QuantumDots.conductance_matrix(::Missing, eig; kwargs...) = missing
function QuantumDots.conductance_matrix(t::Transport, eig; basis=c)
    leads = get_leads(basis, t.parameters...)
    system = t.type(QuantumDots.diagonalize(QuantumDots.OpenSystem(eig, leads)))
    # println(norm(solve(StationaryStateProblem(system))))
    QuantumDots.conductance_matrix(system, .01)
end
function get_leads(c, T, μ, Γ=1)
    l = spatial_labels(c)
    left = QuantumDots.CombinedLead((Γ * c[l[1], :↑]', Γ * c[l[1], :↓]'); T, μ=μ[1])
    right = QuantumDots.CombinedLead((Γ * c[l[end], :↑]', Γ * c[l[end], :↓]'); T, μ=μ[2])
    return (; left, right)
end

