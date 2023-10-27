## Analytic results
ϵodd(t, Ez, θ) = -t + Ez * cos(θ / 2)
ϵodd(; params...) = ϵodd(params[:t], params[:Ez], params[:θ])
ϵeven(Δ, ϕm) = sqrt(2)Δ * cos(ϕm / 2)
ϵeven(; params...) = ϵeven(params[:Δ], params[:ϕm])

Uinf_inds = begin
    cells = [cell_labels(k, c) for k in spatial_labels(c)]
    cellsiteindices = [map(l -> QuantumDots.siteindex(l, c), labels) for labels in cells]
    fs = [f for f in 0:2^length(c)-1 if all(([sum(ind -> QuantumDots._bit(f, ind), siteindices) for siteindices in cellsiteindices]) .< 2)]
    even_fs = [f for f in fs if QuantumDots.parity(f) == 1]
    odd_fs = [f for f in fs if QuantumDots.parity(f) == -1]
    even_inds = map(f -> QuantumDots.focktoind(f, c), even_fs)
    odd_inds = map(f -> QuantumDots.focktoind(f, c), odd_fs)
    return odd_inds, even_inds
end


focknbr(ss::Vararg{AbstractString}) = mapreduce((k, s) -> focknbr(k, s), +, 1:length(ss), ss)
focknbr(s::String) = focknbr(split(s, ",")...)
function focknbr(k::Integer, s::AbstractString)
    s == "0" && return 0
    return QuantumDots.focknbr(QuantumDots.siteindices([(k, Symbol(c)) for c in s], c))
end
index(s) = QuantumDots.focktoind(focknbr(s), c)

qd_state_labels = ["0", "↑", "↓", "↑↓"]
state_labels = [s1 * "," * s2 for (s1, s2) in Base.product(qd_state_labels, qd_state_labels)]
label_to_ind_dict = Dict(map(l -> l => QuantumDots.focktoind(focknbr(l), c), state_labels))

Ψeven(; params...) = Ψeven(params[:θ], params[:ϕp])
function Ψeven(θ, ϕp)
    v = zeros(typeof(complex(θ * ϕp)), 2^length(c))
    v[index("0,0")] = sqrt(2)exp(1im * ϕp / 2)
    v[index("↑,↓")] = cos(θ / 2)
    v[index("↓,↑")] = -cos(θ / 2)
    v[index("↑,↑")] = sin(θ / 2)
    v[index("↓,↓")] = sin(θ / 2)
    v .= v ./ 2
end
Ψodd(; params...) = Ψodd(params[:θ])
function Ψodd(θ)
    v = zeros(typeof(complex(θ)), 2^length(c))
    v[index("↑,0")] = sin(θ / 4)
    v[index("0,↑")] = sin(θ / 4)
    v[index("↓,0")] = cos(θ / 4)
    v[index("0,↓")] = -cos(θ / 4)
    v .= v ./ sqrt(2)
end
