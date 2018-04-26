module JuGP
using JuMP
import Convex
using MathProgBase
using Base.Meta

include("types.jl")
include("operators.jl")
include("expr.jl")
include("solver.jl")

function GPModel(;solver=nothing, method=:LogSumExp)
    m = Model(solver=GPSolver(method, solver))
    m.solvehook = solvehook
    m.ext[:GP] = GPData()
    return m
end

export GPModel

function setdiscretevalues(x::Variable, values)
    m = x.m
    haskey(m.ext, :GP) || error("setdiscretevalues can only be called on GP models")
    vals = convert(Vector{Float64}, collect(values))
    gp = m.ext[:GP]::GPData
    gp.discretevalues[x.col] = vals
    setdiscretevalues
end

export setdiscretevalues
end
