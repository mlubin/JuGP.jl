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
