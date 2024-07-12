using MathLink
#modified from https://github.com/AmplitudeGravity/usingMathLink/
function math2Expr(expr::MathLink.WExpr)
    if expr.head.name == "Times"
        return Expr(:call, :*, map(math2Expr, expr.args)...)
    elseif expr.head.name == "Plus"
        return Expr(:call, :+, map(math2Expr, expr.args)...)
    elseif expr.head.name == "Cos"
        return Expr(:call, :cos, map(math2Expr, expr.args)...)
    elseif expr.head.name == "Sin"
        return Expr(:call, :sin, map(math2Expr, expr.args)...)
    elseif expr.head.name == "Power"
        return Expr(:call, :^, map(math2Expr, expr.args)...)
    elseif expr.head.name == "Rational"
        return Expr(:call, ://, map(math2Expr, expr.args)...)
    elseif expr.head.name == "Sqrt"
        return Expr(:call, :sqrt, map(math2Expr, expr.args)...)
    elseif expr.head.name == "nc"
        return Expr(:call, :*, map(math2Expr, expr.args)...)
    elseif expr.head.name == "List"
        return Expr(:call, :(Base.vect), map(math2Expr, expr.args)...)
    elseif expr.head.name == "a"
        a = Expr(:ref, :(a), map(math2Expr, expr.args[2:end])...)
        if first(expr.args) == 0
            return Expr(:call, :adjoint, a)
        else
            return a
        end
    elseif expr.head.name in ("μ", "ϕ", "ϵ", "Δ")
        return Expr(:ref, Symbol(expr.head.name), map(math2Expr, expr.args)...)
    else
        return Expr(:call, Symbol(expr.head.name), map(math2Expr, expr.args)...)
    end
end
function math2Expr(symb::MathLink.WSymbol)
    if symb.name == "E"
        return :ℯ
    elseif symb.name == "I"
        return :(1im)
    else
        return Symbol(symb.name)
    end
end
math2Expr(num::Number) = num

##
include("perturbation_mathematica_expressions.jl")
zerothorder_JUexp = math2Expr(zerothorder_Wexp)
firstorder_hopping_JUexp = math2Expr(firstorder_hopping_Wexp) |> clipboard
secondorder_N2_nonint_JUexp = math2Expr(secondorder_N2_nonint_Wexp) |> clipboard
secondorder_N3_nonint_JUexp = math2Expr(secondorder_N3_nonint_Wexp) |> clipboard

math2Expr(subs) |> clipboard
math2Expr(subs_homogeneous_Wexp) |> clipboard