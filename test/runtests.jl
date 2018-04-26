#  Copyright 2018, Miles Lubin, Chris Coey, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using JuGP, JuMP
using Base.Test

include(Pkg.dir("JuMP", "test", "solvers.jl"))

methods = [:LogSumExp, :Conic]
cont_solvers = Dict()
cont_solvers[:LogSumExp] = []
cont_solvers[:Conic] = []

ipt && push!(cont_solvers[:LogSumExp], Ipopt.IpoptSolver(print_level=0))
eco && push!(cont_solvers[:Conic], ECOS.ECOSSolver(verbose=false))

include("operators.jl")
include("models.jl")
include("cvx_examples.jl")
