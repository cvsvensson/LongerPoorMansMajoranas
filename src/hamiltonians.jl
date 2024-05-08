function whamiltonian_2site((c1up, c1dn), (c2up, c2dn); t, V, θϕ1, θϕ2, conjugate)
    ms = hopping_rotated_nh(t, (c1up, c1dn), (c2up, c2dn), θϕ1, θϕ2; conjugate)
    if iszero(V)
        return ms
    else
        return ms + V * ((QuantumDots.numberop(c1up) + QuantumDots.numberop(c1dn)) * (QuantumDots.numberop(c2up) + QuantumDots.numberop(c2dn)))
    end
end
function whamiltonian_1site((cup, cdn); ε, Ez, Δ, U, conjugate)
    (ε - Ez) * QuantumDots.numberop(cup) + (ε + Ez) * QuantumDots.numberop(cdn) +
    pairing(Δ, cup, cdn; conjugate) + U * QuantumDots.coulomb(cup, cdn)
end
pairing(Δ, cup, cdn; conjugate) = conjugate ? (Δ * cup' * cdn' + hc) : 2Δ * cup'cdn'

function hopping_rotated_nh(t, (c1up, c1dn), (c2up, c2dn), angles1, angles2; conjugate)
    Ω = QuantumDots.su2_rotation(angles1)' * QuantumDots.su2_rotation(angles2)
    c1 = @SVector [c1up, c1dn]
    c2 = @SVector [c2up, c2dn]
    conjugate ? (t * c1' * Ω * c2 + hc) : 2t * c1' * Ω * c2
end
function whamiltonian(c; ε, Ez, t, Δ, U, V, θ, conjugate)
    M = length(c)
    cell = QuantumDots.cell
    @assert length(cell(1, c)) == 2 "Each unit cell should have two fermions for this hamiltonian"
    N = div(M, 2)
    gv = QuantumDots.getvalue
    h1s = (whamiltonian_1site(cell(j, c); ε=gv(ε, j, N), Ez=gv(Ez, j, N), Δ=gv(Δ, j, N), U=gv(U, j, N), conjugate) for j in 1:N)
    h2s = (whamiltonian_2site(cell(j, c), cell(mod1(j + 1, N), c); t=gv(t, j, N; size=2), V=gv(V, j, N; size=2), θϕ1=(gv(θ, j, N), 0), θϕ2=(gv(θ, mod1(j + 1, N), N), 0), conjugate) for j in 1:N)
    return sum(h1s) + sum(h2s)
end

kitaev_ham(c, ε, Δ, t) = BdGMatrix(QuantumDots.kitaev_hamiltonian(c; μ=-ε, t, Δ); check=false)