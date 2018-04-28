#  Copyright 2016, Miles Lubin and contributors
#  Copyright 2018, Chris Coey and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
    JuMP.initNLP(m)
    m.solvehook = solvehook
    m.ext[:GP] = GPData()
    return m
end

export GPModel

function setdiscretevalues(x::Variable, values)
    m = x.m
    if !haskey(m.ext, :GP)
        error("setdiscretevalues can only be called on a variable of a GP model")
    end
    vals = convert(Vector{Float64}, collect(values))
    gp = m.ext[:GP]::GPData
    gp.discretevalues[x.col] = vals
end

export setdiscretevalues
end
