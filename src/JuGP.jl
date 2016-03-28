module JuGP

using JuMP
import Convex
using MathProgBase
using Base.Meta

# package code goes here
include("types.jl")
include("operators.jl")
include("expr.jl")
include("solver.jl")

function GPModel(;solver=nothing,method=:LogSumExp)
    m = Model(solver=GPSolver(method,solver))
    @defNLParam(m, foo == 1) # to work around NLP resolve warning
    m.solvehook = solvehook
    m.ext[:GP] = GPData()
    return m
end
export GPModel

function setDiscrete(x::Variable, values)
    m = x.m
    haskey(m.ext,:GP) || error("setDiscrete can only be called on GP models")
    vals = convert(Vector{Float64},collect(values))
    gp = m.ext[:GP]::GPData
    gp.discretevalues[x.col] = vals
    nothing
end
export setDiscrete


end # module
