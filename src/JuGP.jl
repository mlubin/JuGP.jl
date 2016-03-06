module JuGP

using JuMP
using MathProgBase
using Base.Meta

# package code goes here
include("types.jl")
include("expr.jl")
include("solver.jl")

end # module
