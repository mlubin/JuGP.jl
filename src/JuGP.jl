module JuGP

using JuMP
using MathProgBase
using Base.Meta

# package code goes here
include("types.jl")
include("operators.jl")
include("expr.jl")
include("solver.jl")

function GPModel()
    m = Model(solver=GPSolver())
    @defNLParam(m, foo == 1) # to work around NLP resolve warning
    m.solvehook = solvehook
    m.ext[:GP] = GPData()
    return m
end
export GPModel

end # module
