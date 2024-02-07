zerothorder_JUexp = :(((adjoint(a[n]) * a[n]) * (μ[n]^2 + -1 * (Vz * sqrt(Δ^2 + μ[n]^2)) + Δ * sqrt(-1 * μ[n] + sqrt(Δ^2 + μ[n]^2)) * sqrt(μ[n] + sqrt(Δ^2 + μ[n]^2)))) * sqrt(Δ^2 + μ[n]^2)^-1);

firstorder_hopping_JUexp = :((tso * (a[n] * a[1+n]) * (ℯ^((1im) * ϕ[1+n]) * sqrt(μ[n] + sqrt(Δ^2 + μ[n]^2)) * sqrt(-1 * μ[1+n] + sqrt(Δ^2 + μ[1+n]^2)) + ℯ^((1im) * ϕ[n]) * sqrt(-1 * μ[n] + sqrt(Δ^2 + μ[n]^2)) * sqrt(μ[1+n] + sqrt(Δ^2 + μ[1+n]^2))) + -1 * (tso * (adjoint(a[n]) * adjoint(a[1+n])) * (ℯ^((1im) * ϕ[n]) * sqrt(μ[n] + sqrt(Δ^2 + μ[n]^2)) * sqrt(-1 * μ[1+n] + sqrt(Δ^2 + μ[1+n]^2)) + ℯ^((1im) * ϕ[1+n]) * sqrt(-1 * μ[n] + sqrt(Δ^2 + μ[n]^2)) * sqrt(μ[1+n] + sqrt(Δ^2 + μ[1+n]^2)))) + t * (adjoint(a[n]) * a[1+n]) * (ℯ^((1im) * ϕ[1+n]) * sqrt(-1 * μ[n] + sqrt(Δ^2 + μ[n]^2)) * sqrt(-1 * μ[1+n] + sqrt(Δ^2 + μ[1+n]^2)) + -1 * (ℯ^((1im) * ϕ[n]) * sqrt(μ[n] + sqrt(Δ^2 + μ[n]^2)) * sqrt(μ[1+n] + sqrt(Δ^2 + μ[1+n]^2)))) + t * (adjoint(a[1+n]) * a[n]) * (ℯ^((1im) * ϕ[n]) * sqrt(-1 * μ[n] + sqrt(Δ^2 + μ[n]^2)) * sqrt(-1 * μ[1+n] + sqrt(Δ^2 + μ[1+n]^2)) + -1 * (ℯ^((1im) * ϕ[1+n]) * sqrt(μ[n] + sqrt(Δ^2 + μ[n]^2)) * sqrt(μ[1+n] + sqrt(Δ^2 + μ[1+n]^2))))) * (2 * ℯ^(((1im) * 2^-1) * (ϕ[n] + ϕ[1+n])) * (Δ^2 + μ[n]^2)^(1 * 4^-1) * (Δ^2 + μ[1+n]^2)^(1 * 4^-1))^-1);

secondorder_N3_nonint_JUexp = :((ℯ^((1im) * ϕ[3]) * (adjoint(a[1]) * a[1]) * (ℯ^((2 * (1im)) * ϕ[1]) * (t^2 + tso^2) * sqrt(-1 * μ[1] + sqrt(Δ^2 + μ[1]^2)) * sqrt(μ[1] + sqrt(Δ^2 + μ[1]^2)) * sqrt(-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[2] + sqrt(Δ^2 + μ[2]^2)) + ℯ^((2 * (1im)) * ϕ[2]) * (t^2 + tso^2) * sqrt(-1 * μ[1] + sqrt(Δ^2 + μ[1]^2)) * sqrt(μ[1] + sqrt(Δ^2 + μ[1]^2)) * sqrt(-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[2] + sqrt(Δ^2 + μ[2]^2)) + 2 * ℯ^((1im) * (ϕ[1] + ϕ[2])) * (t^2 * (-1 * (μ[1] * μ[2]) + sqrt(Δ^2 + μ[1]^2) * sqrt(Δ^2 + μ[2]^2)) + -1 * (tso^2 * (μ[1] * μ[2] + sqrt(Δ^2 + μ[1]^2) * sqrt(Δ^2 + μ[2]^2))))) * sqrt(Δ^2 + μ[3]^2) + ℯ^((1im) * ϕ[1]) * (adjoint(a[3]) * a[3]) * sqrt(Δ^2 + μ[1]^2) * (ℯ^((2 * (1im)) * ϕ[2]) * (t^2 + tso^2) * sqrt(-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(-1 * μ[3] + sqrt(Δ^2 + μ[3]^2)) * sqrt(μ[3] + sqrt(Δ^2 + μ[3]^2)) + ℯ^((2 * (1im)) * ϕ[3]) * (t^2 + tso^2) * sqrt(-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(-1 * μ[3] + sqrt(Δ^2 + μ[3]^2)) * sqrt(μ[3] + sqrt(Δ^2 + μ[3]^2)) + 2 * ℯ^((1im) * (ϕ[2] + ϕ[3])) * (t^2 * (-1 * (μ[2] * μ[3]) + sqrt(Δ^2 + μ[2]^2) * sqrt(Δ^2 + μ[3]^2)) + -1 * (tso^2 * (μ[2] * μ[3] + sqrt(Δ^2 + μ[2]^2) * sqrt(Δ^2 + μ[3]^2))))) + ℯ^(((1im) * 2^-1) * (ϕ[1] + ϕ[3])) * (Δ^2 + μ[1]^2)^(1 * 4^-1) * (Δ^2 + μ[3]^2)^(1 * 4^-1) * ((adjoint(a[1]) * a[3]) * (ℯ^((1im) * (ϕ[1] + ϕ[3])) * (t^2 + -1 * tso^2) * sqrt(μ[1] + sqrt(Δ^2 + μ[1]^2)) * sqrt(-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(-1 * μ[3] + sqrt(Δ^2 + μ[3]^2)) + ℯ^((1im) * (ϕ[2] + ϕ[3])) * sqrt(-1 * μ[1] + sqrt(Δ^2 + μ[1]^2)) * (tso^2 * (-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) + t^2 * (μ[2] + sqrt(Δ^2 + μ[2]^2))) * sqrt(-1 * μ[3] + sqrt(Δ^2 + μ[3]^2)) + ℯ^((2 * (1im)) * ϕ[2]) * (t^2 + -1 * tso^2) * sqrt(-1 * μ[1] + sqrt(Δ^2 + μ[1]^2)) * sqrt(-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[3] + sqrt(Δ^2 + μ[3]^2)) + ℯ^((1im) * (ϕ[1] + ϕ[2])) * sqrt(μ[1] + sqrt(Δ^2 + μ[1]^2)) * (t^2 * (-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) + tso^2 * (μ[2] + sqrt(Δ^2 + μ[2]^2))) * sqrt(μ[3] + sqrt(Δ^2 + μ[3]^2))) + (adjoint(a[3]) * a[1]) * (ℯ^((2 * (1im)) * ϕ[2]) * (t^2 + -1 * tso^2) * sqrt(μ[1] + sqrt(Δ^2 + μ[1]^2)) * sqrt(-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(-1 * μ[3] + sqrt(Δ^2 + μ[3]^2)) + ℯ^((1im) * (ϕ[1] + ϕ[2])) * sqrt(-1 * μ[1] + sqrt(Δ^2 + μ[1]^2)) * (tso^2 * (-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) + t^2 * (μ[2] + sqrt(Δ^2 + μ[2]^2))) * sqrt(-1 * μ[3] + sqrt(Δ^2 + μ[3]^2)) + ℯ^((1im) * (ϕ[1] + ϕ[3])) * (t^2 + -1 * tso^2) * sqrt(-1 * μ[1] + sqrt(Δ^2 + μ[1]^2)) * sqrt(-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[3] + sqrt(Δ^2 + μ[3]^2)) + ℯ^((1im) * (ϕ[2] + ϕ[3])) * sqrt(μ[1] + sqrt(Δ^2 + μ[1]^2)) * (t^2 * (-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) + tso^2 * (μ[2] + sqrt(Δ^2 + μ[2]^2))) * sqrt(μ[3] + sqrt(Δ^2 + μ[3]^2))) + 2 * t * tso * ((a[1] * a[3]) * (ℯ^((1im) * (ϕ[2] + ϕ[3])) * sqrt(μ[1] + sqrt(Δ^2 + μ[1]^2)) * μ[2] * sqrt(-1 * μ[3] + sqrt(Δ^2 + μ[3]^2)) + -1 * (ℯ^((1im) * (ϕ[1] + ϕ[3])) * sqrt(-1 * μ[1] + sqrt(Δ^2 + μ[1]^2)) * sqrt(-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(-1 * μ[3] + sqrt(Δ^2 + μ[3]^2))) + ℯ^((1im) * (ϕ[1] + ϕ[2])) * sqrt(-1 * μ[1] + sqrt(Δ^2 + μ[1]^2)) * μ[2] * sqrt(μ[3] + sqrt(Δ^2 + μ[3]^2)) + ℯ^((2 * (1im)) * ϕ[2]) * sqrt(μ[1] + sqrt(Δ^2 + μ[1]^2)) * sqrt(-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[3] + sqrt(Δ^2 + μ[3]^2))) + (adjoint(a[1]) * adjoint(a[3])) * (-1 * (ℯ^((1im) * (ϕ[1] + ϕ[2])) * sqrt(μ[1] + sqrt(Δ^2 + μ[1]^2)) * μ[2] * sqrt(-1 * μ[3] + sqrt(Δ^2 + μ[3]^2))) + ℯ^((2 * (1im)) * ϕ[2]) * sqrt(-1 * μ[1] + sqrt(Δ^2 + μ[1]^2)) * sqrt(-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(-1 * μ[3] + sqrt(Δ^2 + μ[3]^2)) + -1 * (ℯ^((1im) * (ϕ[2] + ϕ[3])) * sqrt(-1 * μ[1] + sqrt(Δ^2 + μ[1]^2)) * μ[2] * sqrt(μ[3] + sqrt(Δ^2 + μ[3]^2))) + -1 * (ℯ^((1im) * (ϕ[1] + ϕ[3])) * sqrt(μ[1] + sqrt(Δ^2 + μ[1]^2)) * sqrt(-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[3] + sqrt(Δ^2 + μ[3]^2)))))) + (adjoint(a[2]) * a[2]) * (ℯ^((1im) * (2 * ϕ[1] + ϕ[3])) * (t^2 + tso^2) * sqrt(-1 * μ[1] + sqrt(Δ^2 + μ[1]^2)) * sqrt(μ[1] + sqrt(Δ^2 + μ[1]^2)) * sqrt(-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(Δ^2 + μ[3]^2) + ℯ^((1im) * (2 * ϕ[2] + ϕ[3])) * (t^2 + tso^2) * sqrt(-1 * μ[1] + sqrt(Δ^2 + μ[1]^2)) * sqrt(μ[1] + sqrt(Δ^2 + μ[1]^2)) * sqrt(-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(Δ^2 + μ[3]^2) + ℯ^((1im) * (ϕ[1] + 2 * ϕ[2])) * (t^2 + tso^2) * sqrt(Δ^2 + μ[1]^2) * sqrt(-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(-1 * μ[3] + sqrt(Δ^2 + μ[3]^2)) * sqrt(μ[3] + sqrt(Δ^2 + μ[3]^2)) + ℯ^((1im) * (ϕ[1] + 2 * ϕ[3])) * (t^2 + tso^2) * sqrt(Δ^2 + μ[1]^2) * sqrt(-1 * μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(μ[2] + sqrt(Δ^2 + μ[2]^2)) * sqrt(-1 * μ[3] + sqrt(Δ^2 + μ[3]^2)) * sqrt(μ[3] + sqrt(Δ^2 + μ[3]^2)) + ℯ^((1im) * (ϕ[1] + ϕ[2] + ϕ[3])) * (t^2 * (4 * sqrt(Δ^2 + μ[1]^2) * sqrt(Δ^2 + μ[2]^2) * sqrt(Δ^2 + μ[3]^2) + -1 * (2 * μ[2] * (sqrt(Δ^2 + μ[1]^2) * μ[3] + μ[1] * sqrt(Δ^2 + μ[3]^2)))) + -1 * (2 * tso^2 * (2 * sqrt(Δ^2 + μ[1]^2) * sqrt(Δ^2 + μ[2]^2) * sqrt(Δ^2 + μ[3]^2) + μ[2] * (sqrt(Δ^2 + μ[1]^2) * μ[3] + μ[1] * sqrt(Δ^2 + μ[3]^2))))))) * (4 * ℯ^((1im) * (ϕ[1] + ϕ[2] + ϕ[3])) * Vz * sqrt(Δ^2 + μ[1]^2) * sqrt(Δ^2 + μ[2]^2) * sqrt(Δ^2 + μ[3]^2))^-1)

##
##zeroth order
zerothorder_term = @RuntimeGeneratedFunction(:(function zerothorder(n, a, Δ, μ, Vz)
    $zerothorder_JUexp
end))
function zerothorder_perturbation(a; Δ, μ, Ez)
    N = length(a)
    sum(zerothorder_term(n, a, Δ, μ, Ez) for n in 1:N)
end
a = FermionBdGBasis(1:3)
zerothorder_term(3, a, 1, 1:3, 1)

##first order
firstorder_hopping_ex = :(function firstorder_hopping(n, a, Δ, μ, ϕ, tso, t)
    $firstorder_hopping_JUexp
end)
firstorder_hopping_term = @RuntimeGeneratedFunction(firstorder_hopping_ex)
function firstorder_perturbation(a; Δ, μ, δϕ, tso, t)
    N = length(a)
    ϕ = pushfirst!(cumsum(δϕ), 0)
    sum(firstorder_hopping_term(n, a, Δ, μ, ϕ, tso, t) for n in 1:N-1)
end
a = FermionBdGBasis(1:3)
firstorder_perturbation(a; Δ=1, μ=1:3, δϕ=1:2, tso=1, t=1 / 5)
firstorder_hopping_term(1, a, 1.0, rand(3), rand(3), 0.2, 0.2)

## second order
secondorder_hopping_ex = :(function secondorder_hopping(a, Δ, μ, ϕ, tso, t, Vz)
    $secondorder_N3_nonint_JUexp
end);
secondorder_hopping_term = @RuntimeGeneratedFunction(secondorder_hopping_ex)
function secondorder_perturbation(a; Δ, μ, δϕ, tso, t, Ez)
    ϕ = pushfirst!(cumsum(δϕ), 0)
    secondorder_hopping_term(a, Δ, μ, ϕ, tso, t, Ez)
end
a = FermionBdGBasis(1:3)
secondorder_hopping_term(a, 1, 1:3, 1:3, 1, 1, 2)
secondorder_perturbation(a; Δ=1, μ=1:3, δϕ=1:2, tso=1, t=1 / 5, Ez=3)

# Base.:*(a::QuantumDots.BdGFermion, b::QuantumDots.BdGFermion, c::QuantumDots.BdGFermion, d::QuantumDots.BdGFermion) = 0 * (a * b)

##
function perturbative_hamiltonian(a, M=2; Δ, ε, δϕ, t, θ, Ez)
    tso = t * tan(θ / 2)
    _perturbative_hamiltonian(a, M=2; Δ, μ=ε, δϕ, tso, t, Ez)
end
function _perturbative_hamiltonian(a, M=2; Δ, μ, δϕ, tso, t, Ez)
    if M == 0
        zerothorder_perturbation(a; Δ=Δ, μ=μ, Ez=Ez)
    elseif M == 1
        zerothorder_perturbation(a; Δ=Δ, μ=μ, Ez=Ez) + firstorder_perturbation(a; Δ=Δ, μ=μ, δϕ=δϕ, tso=tso, t=t)
    elseif M == 2
        zerothorder_perturbation(a; Δ=Δ, μ=μ, Ez=Ez) + firstorder_perturbation(a; Δ=Δ, μ=μ, δϕ=δϕ, tso=tso, t=t) + secondorder_perturbation(a; Δ=Δ, μ=μ, δϕ=δϕ, tso=tso, t=t, Ez=Ez)
    else
        throw(ArgumentError("M=$M not implemented"))
    end
end
perturbative_hamiltonian(a, 2; Δ=1, μ=1:3, δϕ=1:2, tso=1 / 100000, t=1 / 10000, Ez=3) - zerothorder_perturbation(a; Δ=1, μ=1:3, Ez=3) |> norm
