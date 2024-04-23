subs_homogeneous_JUexp = :(Base.vect(-1 * ((tso * Δ * cos(δϕ * 2^-1)) * sqrt(Δ^2 + ϵ^2)^-1), (tso * Δ * cos(δϕ * 2^-1)) * sqrt(Δ^2 + ϵ^2)^-1, (tso * ϵ * cos(δϕ * 2^-1)) * sqrt(Δ^2 + ϵ^2)^-1 + (1im) * tso * sin(δϕ * 2^-1), (tso * ϵ * cos(δϕ * 2^-1)) * sqrt(Δ^2 + ϵ^2)^-1 + -1 * ((1im) * tso * sin(δϕ * 2^-1)), -1 * ((t * ϵ * cos(δϕ * 2^-1)) * sqrt(Δ^2 + ϵ^2)^-1) + (1im) * t * sin(δϕ * 2^-1), -1 * ((t * Δ * cos(δϕ * 2^-1)) * sqrt(Δ^2 + ϵ^2)^-1), -1 * ((t * Δ * cos(δϕ * 2^-1)) * sqrt(Δ^2 + ϵ^2)^-1), (t * ϵ * cos(δϕ * 2^-1)) * sqrt(Δ^2 + ϵ^2)^-1 + (1im) * t * sin(δϕ * 2^-1)));

subs_JUexp = :(Base.vect((-1 * 2^-1) * ((tso * (ℯ^((1im) * ϕ[1+n]) * sqrt((Δ[n]^2 + ϵ[n] * (ϵ[n] + sqrt(Δ[n]^2 + ϵ[n]^2))) * (Δ[1+n]^2 + ϵ[1+n] * (ϵ[1+n] + -1 * sqrt(Δ[1+n]^2 + ϵ[1+n]^2)))) + ℯ^((1im) * ϕ[n]) * sqrt((Δ[n]^2 + ϵ[n] * (ϵ[n] + -1 * sqrt(Δ[n]^2 + ϵ[n]^2))) * (Δ[1+n]^2 + ϵ[1+n] * (ϵ[1+n] + sqrt(Δ[1+n]^2 + ϵ[1+n]^2)))))) * sqrt((Δ[n]^2 + ϵ[n]^2) * (Δ[1+n]^2 + ϵ[1+n]^2))^-1), (tso * (ℯ^((1im) * ϕ[n]) * sqrt((Δ[n]^2 + ϵ[n] * (ϵ[n] + sqrt(Δ[n]^2 + ϵ[n]^2))) * (Δ[1+n]^2 + ϵ[1+n] * (ϵ[1+n] + -1 * sqrt(Δ[1+n]^2 + ϵ[1+n]^2)))) + ℯ^((1im) * ϕ[1+n]) * sqrt((Δ[n]^2 + ϵ[n] * (ϵ[n] + -1 * sqrt(Δ[n]^2 + ϵ[n]^2))) * (Δ[1+n]^2 + ϵ[1+n] * (ϵ[1+n] + sqrt(Δ[1+n]^2 + ϵ[1+n]^2)))))) * (2 * sqrt((Δ[n]^2 + ϵ[n]^2) * (Δ[1+n]^2 + ϵ[1+n]^2)))^-1, (-1 * 2^-1) * ((tso * (ℯ^((1im) * ϕ[n]) * sqrt((Δ[n]^2 + ϵ[n] * (ϵ[n] + -1 * sqrt(Δ[n]^2 + ϵ[n]^2))) * (Δ[1+n]^2 + ϵ[1+n] * (ϵ[1+n] + -1 * sqrt(Δ[1+n]^2 + ϵ[1+n]^2)))) + -1 * (ℯ^((1im) * ϕ[1+n]) * sqrt((Δ[n]^2 + ϵ[n] * (ϵ[n] + sqrt(Δ[n]^2 + ϵ[n]^2))) * (Δ[1+n]^2 + ϵ[1+n] * (ϵ[1+n] + sqrt(Δ[1+n]^2 + ϵ[1+n]^2))))))) * sqrt((Δ[n]^2 + ϵ[n]^2) * (Δ[1+n]^2 + ϵ[1+n]^2))^-1), (-1 * 2^-1) * ((tso * (ℯ^((1im) * ϕ[1+n]) * sqrt((Δ[n]^2 + ϵ[n] * (ϵ[n] + -1 * sqrt(Δ[n]^2 + ϵ[n]^2))) * (Δ[1+n]^2 + ϵ[1+n] * (ϵ[1+n] + -1 * sqrt(Δ[1+n]^2 + ϵ[1+n]^2)))) + -1 * (ℯ^((1im) * ϕ[n]) * sqrt((Δ[n]^2 + ϵ[n] * (ϵ[n] + sqrt(Δ[n]^2 + ϵ[n]^2))) * (Δ[1+n]^2 + ϵ[1+n] * (ϵ[1+n] + sqrt(Δ[1+n]^2 + ϵ[1+n]^2))))))) * sqrt((Δ[n]^2 + ϵ[n]^2) * (Δ[1+n]^2 + ϵ[1+n]^2))^-1), (t * (ℯ^((1im) * ϕ[1+n]) * sqrt((Δ[n]^2 + ϵ[n] * (ϵ[n] + -1 * sqrt(Δ[n]^2 + ϵ[n]^2))) * (Δ[1+n]^2 + ϵ[1+n] * (ϵ[1+n] + -1 * sqrt(Δ[1+n]^2 + ϵ[1+n]^2)))) + -1 * (ℯ^((1im) * ϕ[n]) * sqrt((Δ[n]^2 + ϵ[n] * (ϵ[n] + sqrt(Δ[n]^2 + ϵ[n]^2))) * (Δ[1+n]^2 + ϵ[1+n] * (ϵ[1+n] + sqrt(Δ[1+n]^2 + ϵ[1+n]^2))))))) * (2 * ℯ^((1im) * ϕ[n]) * sqrt((Δ[n]^2 + ϵ[n]^2) * (Δ[1+n]^2 + ϵ[1+n]^2)))^-1, (-1 * 2^-1) * ((t * (ℯ^((1im) * ϕ[n]) * sqrt((Δ[n]^2 + ϵ[n] * (ϵ[n] + sqrt(Δ[n]^2 + ϵ[n]^2))) * (Δ[1+n]^2 + ϵ[1+n] * (ϵ[1+n] + -1 * sqrt(Δ[1+n]^2 + ϵ[1+n]^2)))) + ℯ^((1im) * ϕ[1+n]) * sqrt((Δ[n]^2 + ϵ[n] * (ϵ[n] + -1 * sqrt(Δ[n]^2 + ϵ[n]^2))) * (Δ[1+n]^2 + ϵ[1+n] * (ϵ[1+n] + sqrt(Δ[1+n]^2 + ϵ[1+n]^2)))))) * (ℯ^((1im) * ϕ[n]) * sqrt((Δ[n]^2 + ϵ[n]^2) * (Δ[1+n]^2 + ϵ[1+n]^2)))^-1), (-1 * 2^-1) * ((t * (ℯ^((1im) * ϕ[1+n]) * sqrt((Δ[n]^2 + ϵ[n] * (ϵ[n] + sqrt(Δ[n]^2 + ϵ[n]^2))) * (Δ[1+n]^2 + ϵ[1+n] * (ϵ[1+n] + -1 * sqrt(Δ[1+n]^2 + ϵ[1+n]^2)))) + ℯ^((1im) * ϕ[n]) * sqrt((Δ[n]^2 + ϵ[n] * (ϵ[n] + -1 * sqrt(Δ[n]^2 + ϵ[n]^2))) * (Δ[1+n]^2 + ϵ[1+n] * (ϵ[1+n] + sqrt(Δ[1+n]^2 + ϵ[1+n]^2)))))) * (ℯ^((1im) * ϕ[n]) * sqrt((Δ[n]^2 + ϵ[n]^2) * (Δ[1+n]^2 + ϵ[1+n]^2)))^-1), (-1 * 2^-1) * ((t * (ℯ^((1im) * ϕ[n]) * sqrt((Δ[n]^2 + ϵ[n] * (ϵ[n] + -1 * sqrt(Δ[n]^2 + ϵ[n]^2))) * (Δ[1+n]^2 + ϵ[1+n] * (ϵ[1+n] + -1 * sqrt(Δ[1+n]^2 + ϵ[1+n]^2)))) + -1 * (ℯ^((1im) * ϕ[1+n]) * sqrt((Δ[n]^2 + ϵ[n] * (ϵ[n] + sqrt(Δ[n]^2 + ϵ[n]^2))) * (Δ[1+n]^2 + ϵ[1+n] * (ϵ[1+n] + sqrt(Δ[1+n]^2 + ϵ[1+n]^2))))))) * (ℯ^((1im) * ϕ[n]) * sqrt((Δ[n]^2 + ϵ[n]^2) * (Δ[1+n]^2 + ϵ[1+n]^2)))^-1)));

zerothorder_JUexp = :((Vz + -1 * sqrt(Δ^2 + ϵ[n]^2)) * (adjoint(a[n]) * a[n]))

firstorder_hopping_JUexp = :((ℯ^((-1 * (1im)) * ϕ[n] + (1im) * ϕ[1+n]) * t * (adjoint(a[n]) * a[1+n]) * sqrt((Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2)) * sqrt((Δ^2 + ϵ[n]^2 + -1 * (ϵ[n] * sqrt(Δ^2 + ϵ[n]^2))) * (Δ^2 + ϵ[1+n]^2 + -1 * (ϵ[1+n] * sqrt(Δ^2 + ϵ[1+n]^2))))) * (2 * (Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2))^-1 + (ℯ^((1im) * ϕ[n] + -1 * ((1im) * ϕ[1+n])) * t * (adjoint(a[1+n]) * a[n]) * sqrt((Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2)) * sqrt((Δ^2 + ϵ[n]^2 + -1 * (ϵ[n] * sqrt(Δ^2 + ϵ[n]^2))) * (Δ^2 + ϵ[1+n]^2 + -1 * (ϵ[1+n] * sqrt(Δ^2 + ϵ[1+n]^2))))) * (2 * (Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2))^-1 + (tso * (adjoint(a[n]) * adjoint(a[1+n])) * sqrt((Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2)) * sqrt((Δ^2 + ϵ[n] * (ϵ[n] + sqrt(Δ^2 + ϵ[n]^2))) * (Δ^2 + ϵ[1+n]^2 + -1 * (ϵ[1+n] * sqrt(Δ^2 + ϵ[1+n]^2))))) * (2 * ℯ^((1im) * ϕ[1+n]) * (Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2))^-1 + -1 * ((ℯ^((1im) * ϕ[1+n]) * tso * (a[n] * a[1+n]) * sqrt((Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2)) * sqrt((Δ^2 + ϵ[n] * (ϵ[n] + sqrt(Δ^2 + ϵ[n]^2))) * (Δ^2 + ϵ[1+n]^2 + -1 * (ϵ[1+n] * sqrt(Δ^2 + ϵ[1+n]^2))))) * (2 * (Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2))^-1) + (tso * (adjoint(a[n]) * adjoint(a[1+n])) * sqrt((Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2)) * sqrt((Δ^2 + ϵ[n]^2 + -1 * (ϵ[n] * sqrt(Δ^2 + ϵ[n]^2))) * (Δ^2 + ϵ[1+n] * (ϵ[1+n] + sqrt(Δ^2 + ϵ[1+n]^2))))) * (2 * ℯ^((1im) * ϕ[n]) * (Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2))^-1 + -1 * ((ℯ^((1im) * ϕ[n]) * tso * (a[n] * a[1+n]) * sqrt((Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2)) * sqrt((Δ^2 + ϵ[n]^2 + -1 * (ϵ[n] * sqrt(Δ^2 + ϵ[n]^2))) * (Δ^2 + ϵ[1+n] * (ϵ[1+n] + sqrt(Δ^2 + ϵ[1+n]^2))))) * (2 * (Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2))^-1) + -1 * ((t * (adjoint(a[n]) * a[1+n]) * sqrt((Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2)) * sqrt((Δ^2 + ϵ[n] * (ϵ[n] + sqrt(Δ^2 + ϵ[n]^2))) * (Δ^2 + ϵ[1+n] * (ϵ[1+n] + sqrt(Δ^2 + ϵ[1+n]^2))))) * (2 * (Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2))^-1) + -1 * ((t * (adjoint(a[1+n]) * a[n]) * sqrt((Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2)) * sqrt((Δ^2 + ϵ[n] * (ϵ[n] + sqrt(Δ^2 + ϵ[n]^2))) * (Δ^2 + ϵ[1+n] * (ϵ[1+n] + sqrt(Δ^2 + ϵ[1+n]^2))))) * (2 * (Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2))^-1))

secondorder_N2_nonint_JUexp = :((ℯ^((-1 * (1im)) * ϕ[n] + (1im) * ϕ[1+n]) * t * (adjoint(a[n]) * a[1+n]) * sqrt((Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2)) * sqrt((Δ^2 + ϵ[n]^2 + -1 * (ϵ[n] * sqrt(Δ^2 + ϵ[n]^2))) * (Δ^2 + ϵ[1+n]^2 + -1 * (ϵ[1+n] * sqrt(Δ^2 + ϵ[1+n]^2))))) * (2 * (Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2))^-1 + (ℯ^((1im) * ϕ[n] + -1 * ((1im) * ϕ[1+n])) * t * (adjoint(a[1+n]) * a[n]) * sqrt((Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2)) * sqrt((Δ^2 + ϵ[n]^2 + -1 * (ϵ[n] * sqrt(Δ^2 + ϵ[n]^2))) * (Δ^2 + ϵ[1+n]^2 + -1 * (ϵ[1+n] * sqrt(Δ^2 + ϵ[1+n]^2))))) * (2 * (Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2))^-1 + (tso * (adjoint(a[n]) * adjoint(a[1+n])) * sqrt((Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2)) * sqrt((Δ^2 + ϵ[n] * (ϵ[n] + sqrt(Δ^2 + ϵ[n]^2))) * (Δ^2 + ϵ[1+n]^2 + -1 * (ϵ[1+n] * sqrt(Δ^2 + ϵ[1+n]^2))))) * (2 * ℯ^((1im) * ϕ[1+n]) * (Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2))^-1 + -1 * ((ℯ^((1im) * ϕ[1+n]) * tso * (a[n] * a[1+n]) * sqrt((Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2)) * sqrt((Δ^2 + ϵ[n] * (ϵ[n] + sqrt(Δ^2 + ϵ[n]^2))) * (Δ^2 + ϵ[1+n]^2 + -1 * (ϵ[1+n] * sqrt(Δ^2 + ϵ[1+n]^2))))) * (2 * (Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2))^-1) + (tso * (adjoint(a[n]) * adjoint(a[1+n])) * sqrt((Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2)) * sqrt((Δ^2 + ϵ[n]^2 + -1 * (ϵ[n] * sqrt(Δ^2 + ϵ[n]^2))) * (Δ^2 + ϵ[1+n] * (ϵ[1+n] + sqrt(Δ^2 + ϵ[1+n]^2))))) * (2 * ℯ^((1im) * ϕ[n]) * (Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2))^-1 + -1 * ((ℯ^((1im) * ϕ[n]) * tso * (a[n] * a[1+n]) * sqrt((Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2)) * sqrt((Δ^2 + ϵ[n]^2 + -1 * (ϵ[n] * sqrt(Δ^2 + ϵ[n]^2))) * (Δ^2 + ϵ[1+n] * (ϵ[1+n] + sqrt(Δ^2 + ϵ[1+n]^2))))) * (2 * (Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2))^-1) + -1 * ((t * (adjoint(a[n]) * a[1+n]) * sqrt((Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2)) * sqrt((Δ^2 + ϵ[n] * (ϵ[n] + sqrt(Δ^2 + ϵ[n]^2))) * (Δ^2 + ϵ[1+n] * (ϵ[1+n] + sqrt(Δ^2 + ϵ[1+n]^2))))) * (2 * (Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2))^-1) + -1 * ((t * (adjoint(a[1+n]) * a[n]) * sqrt((Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2)) * sqrt((Δ^2 + ϵ[n] * (ϵ[n] + sqrt(Δ^2 + ϵ[n]^2))) * (Δ^2 + ϵ[1+n] * (ϵ[1+n] + sqrt(Δ^2 + ϵ[1+n]^2))))) * (2 * (Δ^2 + ϵ[n]^2) * (Δ^2 + ϵ[1+n]^2))^-1))

secondorder_N3_nonint_JUexp = :(
    (
        (-2 * t^2 * (adjoint(a[2]) * a[2])) * (Vz + sqrt(Δ^2 + ϵ[1]^2))^-1 + (2 * tso^2 * (adjoint(a[2]) * a[2])) * (Vz + sqrt(Δ^2 + ϵ[1]^2))^-1 + -1 * ((t^2 * Δ^2 * (adjoint(a[2]) * a[2])) * (ℯ^((1im) * (ϕ[1] + -1 * ϕ[2])) * (Vz + sqrt(Δ^2 + ϵ[1]^2)) * sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[2]^2)))^-1) + -1 * ((ℯ^((1im) * (ϕ[1] + -1 * ϕ[2])) * t^2 * Δ^2 * (adjoint(a[2]) * a[2])) * ((Vz + sqrt(Δ^2 + ϵ[1]^2)) * sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[2]^2)))^-1) + -1 * ((tso^2 * Δ^2 * (adjoint(a[2]) * a[2])) * (ℯ^((1im) * (ϕ[1] + -1 * ϕ[2])) * (Vz + sqrt(Δ^2 + ϵ[1]^2)) * sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[2]^2)))^-1) + -1 * ((ℯ^((1im) * (ϕ[1] + -1 * ϕ[2])) * tso^2 * Δ^2 * (adjoint(a[2]) * a[2])) * ((Vz + sqrt(Δ^2 + ϵ[1]^2)) * sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[2]^2)))^-1) + (2 * t^2 * (adjoint(a[2]) * a[2]) * ϵ[1] * ϵ[2]) * ((Vz + sqrt(Δ^2 + ϵ[1]^2)) * sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[2]^2)))^-1 + (2 * tso^2 * (adjoint(a[2]) * a[2]) * ϵ[1] * ϵ[2]) * ((Vz + sqrt(Δ^2 + ϵ[1]^2)) * sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[2]^2)))^-1 + -1 * ((2 * t^2 * (adjoint(a[1]) * a[1])) * (Vz + sqrt(Δ^2 + ϵ[2]^2))^-1) + (2 * tso^2 * (adjoint(a[1]) * a[1])) * (Vz + sqrt(Δ^2 + ϵ[2]^2))^-1 + -1 * ((2 * t^2 * (adjoint(a[3]) * a[3])) * (Vz + sqrt(Δ^2 + ϵ[2]^2))^-1) + (2 * tso^2 * (adjoint(a[3]) * a[3])) * (Vz + sqrt(Δ^2 + ϵ[2]^2))^-1 + -1 * ((t^2 * Δ^2 * (adjoint(a[1]) * a[1])) * (ℯ^((1im) * (ϕ[1] + -1 * ϕ[2])) * sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[2]^2)) * (Vz + sqrt(Δ^2 + ϵ[2]^2)))^-1) + -1 * ((ℯ^((1im) * (ϕ[1] + -1 * ϕ[2])) * t^2 * Δ^2 * (adjoint(a[1]) * a[1])) * (sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[2]^2)) * (Vz + sqrt(Δ^2 + ϵ[2]^2)))^-1) + -1 * ((tso^2 * Δ^2 * (adjoint(a[1]) * a[1])) * (ℯ^((1im) * (ϕ[1] + -1 * ϕ[2])) * sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[2]^2)) * (Vz + sqrt(Δ^2 + ϵ[2]^2)))^-1) + -1 * ((ℯ^((1im) * (ϕ[1] + -1 * ϕ[2])) * tso^2 * Δ^2 * (adjoint(a[1]) * a[1])) * (sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[2]^2)) * (Vz + sqrt(Δ^2 + ϵ[2]^2)))^-1) + (2 * t^2 * (adjoint(a[1]) * a[1]) * ϵ[1] * ϵ[2]) * (sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[2]^2)) * (Vz + sqrt(Δ^2 + ϵ[2]^2)))^-1 + (2 * tso^2 * (adjoint(a[1]) * a[1]) * ϵ[1] * ϵ[2]) * (sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[2]^2)) * (Vz + sqrt(Δ^2 + ϵ[2]^2)))^-1 + -1 * ((t^2 * Δ^2 * (adjoint(a[3]) * a[3])) * (ℯ^((1im) * (ϕ[2] + -1 * ϕ[3])) * (Vz + sqrt(Δ^2 + ϵ[2]^2)) * sqrt((Δ^2 + ϵ[2]^2) * (Δ^2 + ϵ[3]^2)))^-1) + -1 * ((ℯ^((1im) * (ϕ[2] + -1 * ϕ[3])) * t^2 * Δ^2 * (adjoint(a[3]) * a[3])) * ((Vz + sqrt(Δ^2 + ϵ[2]^2)) * sqrt((Δ^2 + ϵ[2]^2) * (Δ^2 + ϵ[3]^2)))^-1) + -1 * ((tso^2 * Δ^2 * (adjoint(a[3]) * a[3])) * (ℯ^((1im) * (ϕ[2] + -1 * ϕ[3])) * (Vz + sqrt(Δ^2 + ϵ[2]^2)) * sqrt((Δ^2 + ϵ[2]^2) * (Δ^2 + ϵ[3]^2)))^-1) + -1 * ((ℯ^((1im) * (ϕ[2] + -1 * ϕ[3])) * tso^2 * Δ^2 * (adjoint(a[3]) * a[3])) * ((Vz + sqrt(Δ^2 + ϵ[2]^2)) * sqrt((Δ^2 + ϵ[2]^2) * (Δ^2 + ϵ[3]^2)))^-1) + (2 * t^2 * (adjoint(a[3]) * a[3]) * ϵ[2] * ϵ[3]) * ((Vz + sqrt(Δ^2 + ϵ[2]^2)) * sqrt((Δ^2 + ϵ[2]^2) * (Δ^2 + ϵ[3]^2)))^-1 + (2 * tso^2 * (adjoint(a[3]) * a[3]) * ϵ[2] * ϵ[3]) * ((Vz + sqrt(Δ^2 + ϵ[2]^2)) * sqrt((Δ^2 + ϵ[2]^2) * (Δ^2 + ϵ[3]^2)))^-1 + -1 * ((tso^2 * (adjoint(a[1]) * a[3]) * (1 + -1 * (ϵ[2] * sqrt(Δ^2 + ϵ[2]^2)^-1)) * sqrt((1 + -1 * (ϵ[1] * sqrt(Δ^2 + ϵ[1]^2)^-1)) * (1 + -1 * (ϵ[3] * sqrt(Δ^2 + ϵ[3]^2)^-1)))) * (ℯ^((1im) * (ϕ[1] + -1 * ϕ[3])) * (Vz + sqrt(Δ^2 + ϵ[2]^2)))^-1) + -1 * ((ℯ^((1im) * (ϕ[1] + -1 * ϕ[3])) * tso^2 * (adjoint(a[3]) * a[1]) * (1 + -1 * (ϵ[2] * sqrt(Δ^2 + ϵ[2]^2)^-1)) * sqrt((1 + -1 * (ϵ[1] * sqrt(Δ^2 + ϵ[1]^2)^-1)) * (1 + -1 * (ϵ[3] * sqrt(Δ^2 + ϵ[3]^2)^-1)))) * (Vz + sqrt(Δ^2 + ϵ[2]^2))^-1) + (t * tso * (adjoint(a[1]) * adjoint(a[3])) * (1 + -1 * (ϵ[2] * sqrt(Δ^2 + ϵ[2]^2)^-1)) * sqrt((ϵ[1] + sqrt(Δ^2 + ϵ[1]^2)) * (1 + -1 * (ϵ[3] * sqrt(Δ^2 + ϵ[3]^2)^-1)))) * (ℯ^((1im) * ϕ[3]) * (Δ^2 + ϵ[1]^2)^(1 * 4^-1) * (Vz + sqrt(Δ^2 + ϵ[2]^2)))^-1 + -1 * ((ℯ^((1im) * ϕ[3]) * t * tso * (a[1] * a[3]) * (1 + -1 * (ϵ[2] * sqrt(Δ^2 + ϵ[2]^2)^-1)) * sqrt((ϵ[1] + sqrt(Δ^2 + ϵ[1]^2)) * (1 + -1 * (ϵ[3] * sqrt(Δ^2 + ϵ[3]^2)^-1)))) * ((Δ^2 + ϵ[1]^2)^(1 * 4^-1) * (Vz + sqrt(Δ^2 + ϵ[2]^2)))^-1) + (2 * t * tso * Δ * (adjoint(a[1]) * adjoint(a[3])) * sqrt(((1 + -1 * (ϵ[1] * sqrt(Δ^2 + ϵ[1]^2)^-1)) * (1 + -1 * (ϵ[3] * sqrt(Δ^2 + ϵ[3]^2)^-1))) * (Δ^2 + ϵ[2]^2)^-1)) * (ℯ^((1im) * (ϕ[1] + -1 * ϕ[2] + ϕ[3])) * (Vz + sqrt(Δ^2 + ϵ[2]^2)))^-1 + -1 * ((2 * ℯ^((1im) * (ϕ[1] + -1 * ϕ[2] + ϕ[3])) * t * tso * Δ * (a[1] * a[3]) * sqrt(((1 + -1 * (ϵ[1] * sqrt(Δ^2 + ϵ[1]^2)^-1)) * (1 + -1 * (ϵ[3] * sqrt(Δ^2 + ϵ[3]^2)^-1))) * (Δ^2 + ϵ[2]^2)^-1)) * (Vz + sqrt(Δ^2 + ϵ[2]^2))^-1) + -1 * ((t^2 * (adjoint(a[1]) * a[3]) * (ϵ[2] + sqrt(Δ^2 + ϵ[2]^2)) * sqrt(((1 + -1 * (ϵ[1] * sqrt(Δ^2 + ϵ[1]^2)^-1)) * (1 + -1 * (ϵ[3] * sqrt(Δ^2 + ϵ[3]^2)^-1))) * (Δ^2 + ϵ[2]^2)^-1)) * (ℯ^((1im) * (ϕ[1] + -1 * ϕ[3])) * (Vz + sqrt(Δ^2 + ϵ[2]^2)))^-1) + -1 * ((ℯ^((1im) * (ϕ[1] + -1 * ϕ[3])) * t^2 * (adjoint(a[3]) * a[1]) * (ϵ[2] + sqrt(Δ^2 + ϵ[2]^2)) * sqrt(((1 + -1 * (ϵ[1] * sqrt(Δ^2 + ϵ[1]^2)^-1)) * (1 + -1 * (ϵ[3] * sqrt(Δ^2 + ϵ[3]^2)^-1))) * (Δ^2 + ϵ[2]^2)^-1)) * (Vz + sqrt(Δ^2 + ϵ[2]^2))^-1) + -1 * ((2 * t^2 * (adjoint(a[2]) * a[2])) * (Vz + sqrt(Δ^2 + ϵ[3]^2))^-1) + (2 * tso^2 * (adjoint(a[2]) * a[2])) * (Vz + sqrt(Δ^2 + ϵ[3]^2))^-1 + -1 * ((t^2 * Δ^2 * (adjoint(a[2]) * a[2])) * (ℯ^((1im) * (ϕ[2] + -1 * ϕ[3])) * sqrt((Δ^2 + ϵ[2]^2) * (Δ^2 + ϵ[3]^2)) * (Vz + sqrt(Δ^2 + ϵ[3]^2)))^-1) + -1 * ((ℯ^((1im) * (ϕ[2] + -1 * ϕ[3])) * t^2 * Δ^2 * (adjoint(a[2]) * a[2])) * (sqrt((Δ^2 + ϵ[2]^2) * (Δ^2 + ϵ[3]^2)) * (Vz + sqrt(Δ^2 + ϵ[3]^2)))^-1) + -1 * ((tso^2 * Δ^2 * (adjoint(a[2]) * a[2])) * (ℯ^((1im) * (ϕ[2] + -1 * ϕ[3])) * sqrt((Δ^2 + ϵ[2]^2) * (Δ^2 + ϵ[3]^2)) * (Vz + sqrt(Δ^2 + ϵ[3]^2)))^-1) + -1 * ((ℯ^((1im) * (ϕ[2] + -1 * ϕ[3])) * tso^2 * Δ^2 * (adjoint(a[2]) * a[2])) * (sqrt((Δ^2 + ϵ[2]^2) * (Δ^2 + ϵ[3]^2)) * (Vz + sqrt(Δ^2 + ϵ[3]^2)))^-1) + (2 * t^2 * (adjoint(a[2]) * a[2]) * ϵ[2] * ϵ[3]) * (sqrt((Δ^2 + ϵ[2]^2) * (Δ^2 + ϵ[3]^2)) * (Vz + sqrt(Δ^2 + ϵ[3]^2)))^-1 + (2 * tso^2 * (adjoint(a[2]) * a[2]) * ϵ[2] * ϵ[3]) * (sqrt((Δ^2 + ϵ[2]^2) * (Δ^2 + ϵ[3]^2)) * (Vz + sqrt(Δ^2 + ϵ[3]^2)))^-1 + -1 * ((t * tso * (adjoint(a[1]) * adjoint(a[3])) * (ϵ[2] + sqrt(Δ^2 + ϵ[2]^2)) * sqrt(((ϵ[1] + sqrt(Δ^2 + ϵ[1]^2)) * (-1 * ϵ[3] + sqrt(Δ^2 + ϵ[3]^2))) * (Δ^2 + ϵ[2]^2)^-1)) * (ℯ^((1im) * ϕ[3]) * (Vz + sqrt(Δ^2 + ϵ[2]^2)) * ((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[3]^2))^(1 * 4^-1))^-1) + (ℯ^((1im) * ϕ[3]) * t * tso * (a[1] * a[3]) * (ϵ[2] + sqrt(Δ^2 + ϵ[2]^2)) * sqrt(((ϵ[1] + sqrt(Δ^2 + ϵ[1]^2)) * (-1 * ϵ[3] + sqrt(Δ^2 + ϵ[3]^2))) * (Δ^2 + ϵ[2]^2)^-1)) * ((Vz + sqrt(Δ^2 + ϵ[2]^2)) * ((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[3]^2))^(1 * 4^-1))^-1 + (t * tso * (adjoint(a[1]) * adjoint(a[3])) * (1 + -1 * (ϵ[2] * sqrt(Δ^2 + ϵ[2]^2)^-1)) * sqrt((1 + -1 * (ϵ[1] * sqrt(Δ^2 + ϵ[1]^2)^-1)) * (ϵ[3] + sqrt(Δ^2 + ϵ[3]^2)))) * (ℯ^((1im) * ϕ[1]) * (Vz + sqrt(Δ^2 + ϵ[2]^2)) * (Δ^2 + ϵ[3]^2)^(1 * 4^-1))^-1 + -1 * ((ℯ^((1im) * ϕ[1]) * t * tso * (a[1] * a[3]) * (1 + -1 * (ϵ[2] * sqrt(Δ^2 + ϵ[2]^2)^-1)) * sqrt((1 + -1 * (ϵ[1] * sqrt(Δ^2 + ϵ[1]^2)^-1)) * (ϵ[3] + sqrt(Δ^2 + ϵ[3]^2)))) * ((Vz + sqrt(Δ^2 + ϵ[2]^2)) * (Δ^2 + ϵ[3]^2)^(1 * 4^-1))^-1) + -1 * ((t^2 * (adjoint(a[1]) * a[3]) * (1 + -1 * (ϵ[2] * sqrt(Δ^2 + ϵ[2]^2)^-1)) * sqrt((ϵ[1] + sqrt(Δ^2 + ϵ[1]^2)) * (ϵ[3] + sqrt(Δ^2 + ϵ[3]^2)))) * ((Vz + sqrt(Δ^2 + ϵ[2]^2)) * ((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[3]^2))^(1 * 4^-1))^-1) + -1 * ((t^2 * (adjoint(a[3]) * a[1]) * (1 + -1 * (ϵ[2] * sqrt(Δ^2 + ϵ[2]^2)^-1)) * sqrt((ϵ[1] + sqrt(Δ^2 + ϵ[1]^2)) * (ϵ[3] + sqrt(Δ^2 + ϵ[3]^2)))) * ((Vz + sqrt(Δ^2 + ϵ[2]^2)) * ((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[3]^2))^(1 * 4^-1))^-1) + -1 * ((t * tso * (adjoint(a[1]) * adjoint(a[3])) * (ϵ[2] + sqrt(Δ^2 + ϵ[2]^2)) * sqrt(((-1 * ϵ[1] + sqrt(Δ^2 + ϵ[1]^2)) * (ϵ[3] + sqrt(Δ^2 + ϵ[3]^2))) * (Δ^2 + ϵ[2]^2)^-1)) * (ℯ^((1im) * ϕ[1]) * (Vz + sqrt(Δ^2 + ϵ[2]^2)) * ((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[3]^2))^(1 * 4^-1))^-1) + (ℯ^((1im) * ϕ[1]) * t * tso * (a[1] * a[3]) * (ϵ[2] + sqrt(Δ^2 + ϵ[2]^2)) * sqrt(((-1 * ϵ[1] + sqrt(Δ^2 + ϵ[1]^2)) * (ϵ[3] + sqrt(Δ^2 + ϵ[3]^2))) * (Δ^2 + ϵ[2]^2)^-1)) * ((Vz + sqrt(Δ^2 + ϵ[2]^2)) * ((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[3]^2))^(1 * 4^-1))^-1 + -1 * ((tso^2 * (adjoint(a[1]) * a[3]) * (ϵ[2] + sqrt(Δ^2 + ϵ[2]^2)) * sqrt(((ϵ[1] + sqrt(Δ^2 + ϵ[1]^2)) * (ϵ[3] + sqrt(Δ^2 + ϵ[3]^2))) * (Δ^2 + ϵ[2]^2)^-1)) * ((Vz + sqrt(Δ^2 + ϵ[2]^2)) * ((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[3]^2))^(1 * 4^-1))^-1) + -1 * ((tso^2 * (adjoint(a[3]) * a[1]) * (ϵ[2] + sqrt(Δ^2 + ϵ[2]^2)) * sqrt(((ϵ[1] + sqrt(Δ^2 + ϵ[1]^2)) * (ϵ[3] + sqrt(Δ^2 + ϵ[3]^2))) * (Δ^2 + ϵ[2]^2)^-1)) * ((Vz + sqrt(Δ^2 + ϵ[2]^2)) * ((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[3]^2))^(1 * 4^-1))^-1) + -1 * ((ℯ^((1im) * (ϕ[2] + -1 * ϕ[3])) * t^2 * Δ * (adjoint(a[3]) * a[1]) * sqrt(1 + ϵ[1] * sqrt(Δ^2 + ϵ[1]^2)^-1 + -1 * (ϵ[3] * sqrt(Δ^2 + ϵ[3]^2)^-1) + -1 * ((ϵ[1] * ϵ[3]) * sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[3]^2))^-1))) * (Δ^2 + ϵ[2]^2 + Vz * sqrt(Δ^2 + ϵ[2]^2))^-1) + (ℯ^((1im) * (ϕ[2] + -1 * ϕ[3])) * tso^2 * Δ * (adjoint(a[3]) * a[1]) * sqrt(1 + ϵ[1] * sqrt(Δ^2 + ϵ[1]^2)^-1 + -1 * (ϵ[3] * sqrt(Δ^2 + ϵ[3]^2)^-1) + -1 * ((ϵ[1] * ϵ[3]) * sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[3]^2))^-1))) * (Δ^2 + ϵ[2]^2 + Vz * sqrt(Δ^2 + ϵ[2]^2))^-1 + -1 * ((t^2 * Δ * (adjoint(a[1]) * a[3]) * sqrt((1 + ϵ[1] * sqrt(Δ^2 + ϵ[1]^2)^-1 + -1 * (ϵ[3] * sqrt(Δ^2 + ϵ[3]^2)^-1) + -1 * ((ϵ[1] * ϵ[3]) * sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[3]^2))^-1)) * (Δ^2 + ϵ[2]^2)^-1)) * (ℯ^((1im) * (ϕ[2] + -1 * ϕ[3])) * (Vz + sqrt(Δ^2 + ϵ[2]^2)))^-1) + (tso^2 * Δ * (adjoint(a[1]) * a[3]) * sqrt((1 + ϵ[1] * sqrt(Δ^2 + ϵ[1]^2)^-1 + -1 * (ϵ[3] * sqrt(Δ^2 + ϵ[3]^2)^-1) + -1 * ((ϵ[1] * ϵ[3]) * sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[3]^2))^-1)) * (Δ^2 + ϵ[2]^2)^-1)) * (ℯ^((1im) * (ϕ[2] + -1 * ϕ[3])) * (Vz + sqrt(Δ^2 + ϵ[2]^2)))^-1 + -1 * ((ℯ^((1im) * (ϕ[1] + -1 * ϕ[2])) * t^2 * Δ * (adjoint(a[3]) * a[1]) * sqrt(1 + -1 * (ϵ[1] * sqrt(Δ^2 + ϵ[1]^2)^-1) + ϵ[3] * sqrt(Δ^2 + ϵ[3]^2)^-1 + -1 * ((ϵ[1] * ϵ[3]) * sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[3]^2))^-1))) * (Δ^2 + ϵ[2]^2 + Vz * sqrt(Δ^2 + ϵ[2]^2))^-1) +
        (ℯ^((1im) * (ϕ[1] + -1 * ϕ[2])) * tso^2 * Δ * (adjoint(a[3]) * a[1]) * sqrt(1 + -1 * (ϵ[1] * sqrt(Δ^2 + ϵ[1]^2)^-1) + ϵ[3] * sqrt(Δ^2 + ϵ[3]^2)^-1 + -1 * ((ϵ[1] * ϵ[3]) * sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[3]^2))^-1))) * (Δ^2 + ϵ[2]^2 + Vz * sqrt(Δ^2 + ϵ[2]^2))^-1 + -1 * ((t^2 * Δ * (adjoint(a[1]) * a[3]) * sqrt((1 + -1 * (ϵ[1] * sqrt(Δ^2 + ϵ[1]^2)^-1) + ϵ[3] * sqrt(Δ^2 + ϵ[3]^2)^-1 + -1 * ((ϵ[1] * ϵ[3]) * sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[3]^2))^-1)) * (Δ^2 + ϵ[2]^2)^-1)) * (ℯ^((1im) * (ϕ[1] + -1 * ϕ[2])) * (Vz + sqrt(Δ^2 + ϵ[2]^2)))^-1) + (tso^2 * Δ * (adjoint(a[1]) * a[3]) * sqrt((1 + -1 * (ϵ[1] * sqrt(Δ^2 + ϵ[1]^2)^-1) + ϵ[3] * sqrt(Δ^2 + ϵ[3]^2)^-1 + -1 * ((ϵ[1] * ϵ[3]) * sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[3]^2))^-1)) * (Δ^2 + ϵ[2]^2)^-1)) * (ℯ^((1im) * (ϕ[1] + -1 * ϕ[2])) * (Vz + sqrt(Δ^2 + ϵ[2]^2)))^-1 + (2 * ℯ^((1im) * ϕ[2]) * t * tso * Δ * (a[1] * a[3]) * sqrt(1 + ϵ[1] * sqrt(Δ^2 + ϵ[1]^2)^-1 + ϵ[3] * sqrt(Δ^2 + ϵ[3]^2)^-1 + (ϵ[1] * ϵ[3]) * sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[3]^2))^-1)) * (Δ^2 + ϵ[2]^2 + Vz * sqrt(Δ^2 + ϵ[2]^2))^-1 + -1 * ((2 * t * tso * Δ * (adjoint(a[1]) * adjoint(a[3])) * sqrt((1 + ϵ[1] * sqrt(Δ^2 + ϵ[1]^2)^-1 + ϵ[3] * sqrt(Δ^2 + ϵ[3]^2)^-1 + (ϵ[1] * ϵ[3]) * sqrt((Δ^2 + ϵ[1]^2) * (Δ^2 + ϵ[3]^2))^-1) * (Δ^2 + ϵ[2]^2)^-1)) * (ℯ^((1im) * ϕ[2]) * (Vz + sqrt(Δ^2 + ϵ[2]^2)))^-1)
    ) * 4^-1
)

##
##zeroth order
zerothorder_term = @RuntimeGeneratedFunction(:(function zerothorder(n, a, Δ, ϵ, Vz)
    $zerothorder_JUexp
end))
function zerothorder_perturbation(a; Δ, ε, Ez)
    N = length(a)
    sum(zerothorder_term(n, a, Δ, ε, Ez) for n in 1:N)
end
a = FermionBdGBasis(1:3)
zerothorder_term(3, a, 1, [1, 2, 3], 1)

##first order
firstorder_hopping_ex = :(function firstorder_hopping(n, a, Δ, ϵ, ϕ, tso, t)
    $firstorder_hopping_JUexp
end)
firstorder_hopping_term = @RuntimeGeneratedFunction(firstorder_hopping_ex)
function firstorder_perturbation(a; Δ, ε, δϕ, tso, t)
    N = length(a)
    ϕ = pushfirst!(cumsum(δϕ), 0)
    sum(firstorder_hopping_term(n, a, Δ, ε, ϕ, tso, t) for n in 1:N-1)
end
a = FermionBdGBasis(1:3)
firstorder_perturbation(a; Δ=1, ε=[1, 2, 3], δϕ=1:2, tso=1, t=1 / 5)
firstorder_hopping_term(1, a, 1.0, rand(3), rand(3), 0.2, 0.2)

## second order
secondorder_hopping_ex = :(function secondorder_hopping(a, Δ, ϵ, ϕ, tso, t, Vz)
    if length(a) == 3
        return $secondorder_N3_nonint_JUexp
    elseif length(a) == 2
        return $secondorder_N2_nonint_JUexp
    else
        throw(ArgumentError("Only N=2 or N=3 supported"))
    end
end);
secondorder_hopping_term = @RuntimeGeneratedFunction(secondorder_hopping_ex)
function secondorder_perturbation(a; Δ, ε, δϕ, tso, t, Ez)
    ϕ = pushfirst!(cumsum(δϕ), 0)
    secondorder_hopping_term(a, Δ, ε, ϕ, tso, t, Ez)
end
a = FermionBdGBasis(1:3)
secondorder_hopping_term(a, 1, 1:3, 1:3, 1, 1, 2)
secondorder_perturbation(a; Δ=1, ε=[1, 2, 3], δϕ=1:2, tso=1, t=1 / 5, Ez=3)

# Base.:*(a::QuantumDots.BdGFermion, b::QuantumDots.BdGFermion, c::QuantumDots.BdGFermion, d::QuantumDots.BdGFermion) = 0 * (a * b)

##
function theta_to_ts(θ::Number, t)
    N = sqrt(1 + tan(θ / 2)^2)
    (t / N, t * tan(θ / 2) / N)
end
function theta_to_ts(θ::QuantumDots.DiffChainParameter, t)
    theta_to_ts(θ.value, t)
end
function perturbative_hamiltonian_old(a, M; Δ, ε, δϕ, t, θ, Ez)
    t, tso = theta_to_ts(θ, t)
    _perturbative_hamiltonian(a, M; Δ, ε, δϕ, tso, t, Ez)
end
function perturbative_hamiltonian_terms(a; Δ, ε, δϕ, t, θ, Ez)
    t, tso = theta_to_ts(θ, t)
    _perturbative_hamiltonian_terms(a; Δ, ε, δϕ, tso, t, Ez)
end
function _perturbative_hamiltonian_terms(a; Δ, ε, δϕ, tso, t, Ez)
    [zerothorder_perturbation(a; Δ, ε, Ez),
        firstorder_perturbation(a; Δ, ε, δϕ, tso, t),
        secondorder_perturbation(a; Δ, ε, δϕ, tso, t, Ez)]
end
function _perturbative_hamiltonian(a, M; Δ, ε, δϕ, tso, t, Ez)
    perts = _perturbative_hamiltonian_terms(a; Δ, ε, δϕ, tso, t, Ez)
    return sum(perts[1:M+1])
end
perturbative_hamiltonian_old(a, 2; Δ=1, ε=[1, 2, 3], δϕ=1:2, t=1 / 100000, θ=2.7, Ez=3) - zerothorder_perturbation(a; Δ=1, ε=[1, 2, 3], Ez=3) |> norm

##subs
subs_term = @RuntimeGeneratedFunction(:(function subs_JUexp_func(n, Δ, ϵ, t, tso, ϕ)
    $subs_JUexp
end))
function perturbative_coeffs(n; Δ, ε, δϕ, θ, t)
    ϕ = pushfirst!(cumsum(δϕ), zero(eltype(δϕ)))
    t, tso = theta_to_ts(θ, t)
    subs_term(n, Δ, ε, t, tso, ϕ)
end
perturbative_coeffs(2; Δ=[1, 1, 1], ε=[1, 2, 3], θ=0.2, δϕ=1:2, t=1)


##
function perturbative_hamiltonian(a, M; Δ, ε, δϕ, t, θ, Ez)
    N = length(a)
    coeffs = stack(perturbative_coeffs(n; Δ, ε, θ, δϕ, t) for n in 1:N-1)
    Δ_aa, Δ_bb, Δ_ab, Δ_ba, t_aa, t_ab, t_ba, t_bb = collect(eachrow(coeffs))
    perts = [zeroth_order_perturbative_hamiltonian(a; Δ, ε, Ez),
        first_order_perturbative_hamiltonian(a; Δ_aa, t_aa),
        second_order_perturbative_hamiltonian(a; Δ_ba, Δ_ab, t_ba, t_ab, Ez, ε, Δ)]
    sum(perts[1:M+1])
end
function zeroth_order_perturbative_hamiltonian(a; Ez, ε, Δ)
    N = length(a)
    sum((Ez - sqrt(Δ[n]^2 + ε[n]^2)) * a[n]'a[n] for n in 1:N)
end
function first_order_perturbative_hamiltonian(a; t_aa, Δ_aa)
    N = length(a)
    sum(a[n]'a[n+1] * t_aa[n] + Δ_aa[n] * a[n] * a[n+1] for n in 1:N-1) + hc
end

function second_order_coeffs(n; Δ_ba, Δ_ab, t_ba, t_ab, Ez, ε, Δ)
    E1 = Ez + sqrt(Δ[n]^2 + ε[n]^2)
    E2 = Ez + sqrt(Δ[n+1]^2 + ε[n+1]^2)
    ε_ba = (Δ_ba[n]' * Δ_ba[n] - t_ba[n]' * t_ba[n]) / E1
    ε_ab = (Δ_ab[n]' * Δ_ab[n] - t_ab[n]' * t_ab[n]) / E2
    t_nn = (-Δ_ab[n]' * Δ_ba[n+1] - t_ab[n] * t_ba[n+1]) / E2
    Δ_nn = (Δ_ba[n+1]' * t_ab[n] + t_ba[n+1]' * Δ_ab[n]') / E2
    (; E1, E2, ε_ba, ε_ab, t_nn, Δ_nn)
end
function second_order_perturbative_hamiltonian(a; Δ_ba, Δ_ab, t_ba, t_ab, Ez, ε, Δ)
    N = length(a)
    E = [Ez + sqrt(Δ[n]^2 + ε[n]^2) for n in 1:N]
    H = 0 * (a[1] * a[1])

    for n in 1:N-1
        ε_ba = (Δ_ba[n]' * Δ_ba[n] - t_ba[n]' * t_ba[n]) / E[n]
        ε_ab = (Δ_ab[n]' * Δ_ab[n] - t_ab[n]' * t_ab[n]) / E[n+1]
        H += ε_ba * a[n+1]'a[n+1] + ε_ab * a[n]'a[n]
    end

    for n in 1:N-2
        t_nn = (-Δ_ab[n]' * Δ_ba[n+1] - t_ab[n] * t_ba[n+1]) / E[n+1]
        Δ_nn = (Δ_ba[n+1]' * t_ab[n] + t_ba[n+1]' * Δ_ab[n]') / E[n+1]
        H += (t_nn * a[n]' * a[n+2] + Δ_nn * a[n]' * a[n+2]' + hc)
    end
    return H
end
perturbative_hamiltonian_old(a, 2; Δ=1, ε=[1, 2, 3], δϕ=1:2, t=0.5, θ=2.7, Ez=3) -
perturbative_hamiltonian(a, 2; Δ=[1, 1, 1], ε=[1, 2, 3], δϕ=1:2, t=0.5, θ=2.7, Ez=3) |> norm

##
##subs
subs_homogeneous = @RuntimeGeneratedFunction(:(function subs_homogeneous_func(Δ, ϵ, t, tso, δϕ)
    $subs_homogeneous_JUexp
end))
function perturbative_coeffs_homogeneous(; Δ, ε, δϕ, θ, t)
    t, tso = theta_to_ts(θ, t)
    subs_homogeneous(Δ, ε, t, tso, δϕ)
end
perturbative_coeffs_homogeneous(; Δ=1, ε=0.1, θ=0.2, δϕ=1.2, t=1)

function perturbative_hamiltonian_homogeneous(a, M; Δ, ε, δϕ, t, θ, Ez)
    Δ_aa, Δ_bb, Δ_ab, Δ_ba, t_aa, t_ab, t_ba, t_bb = perturbative_coeffs_homogeneous(; Δ, ε, θ, δϕ, t)
    perts = [zeroth_order_perturbative_hamiltonian_homogeneous(a; Δ, ε, Ez),
        first_order_perturbative_hamiltonian_homogeneous(a; Δ_aa, t_aa),
        second_order_perturbative_hamiltonian_homogeneous(a; Δ_ba, Δ_ab, t_ba, t_ab, Ez, ε, Δ)]
    sum(perts[1:M+1])
end
function zeroth_order_perturbative_hamiltonian_homogeneous(a; Ez, ε, Δ)
    N = length(a)
    sum((Ez - sqrt(Δ^2 + ε^2)) * a[n]'a[n] for n in 1:N)
end
function first_order_perturbative_hamiltonian_homogeneous(a; t_aa, Δ_aa)
    N = length(a)
    sum(a[n]'a[n+1] * t_aa + Δ_aa * a[n] * a[n+1] for n in 1:N-1) + hc
end

function second_order_coeffs_homogeneous(; Δ_ba, Δ_ab, t_ba, t_ab, Ez, ε, Δ)
    E = Ez + sqrt(Δ^2 + ε^2)
    ε_ba = (Δ_ba' * Δ_ba - t_ba' * t_ba) / E
    ε_ab = (Δ_ab' * Δ_ab - t_ab' * t_ab) / E
    t_nn = (-Δ_ab' * Δ_ba - t_ab * t_ba) / E
    Δ_nn = (Δ_ba' * t_ab + t_ba' * Δ_ab') / E
    (; E, ε_ba, ε_ab, t_nn, Δ_nn)
end
function second_order_perturbative_hamiltonian_homogeneous(a; Δ_ba, Δ_ab, t_ba, t_ab, Ez, ε, Δ)
    N = length(a)
    E = Ez + sqrt(Δ^2 + ε^2)
    (; E, ε_ba, ε_ab, t_nn, Δ_nn) = second_order_coeffs_homogeneous(; Δ_ba, Δ_ab, t_ba, t_ab, Ez, ε, Δ)
    H = 0 * (a[1] * a[1])

    for n in 1:N-1
        H += ε_ba * a[n+1]'a[n+1] + ε_ab * a[n]'a[n]
    end

    for n in 1:N-2
        H += (t_nn * a[n]' * a[n+2] + Δ_nn * a[n]' * a[n+2]' + hc)
    end
    return H
end

norm(perturbative_hamiltonian_old(a, 2; Δ=1, ε=[1, 1, 1], δϕ=[1, 1], t=0.5, θ=2.7, Ez=3)) -
norm(perturbative_hamiltonian_homogeneous(a, 2; Δ=1, ε=1, δϕ=1, t=0.5, θ=2.7, Ez=3))
norm(perturbative_hamiltonian(a, 2; Δ=[1, 1, 1], ε=[1, 1, 1], δϕ=[1, 1], t=0.5, θ=2.7, Ez=3)) -
norm(perturbative_hamiltonian_homogeneous(a, 2; Δ=1, ε=1, δϕ=1, t=0.5, θ=2.7, Ez=3))