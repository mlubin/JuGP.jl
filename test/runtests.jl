#  Copyright 2016, Miles Lubin and contributors
#  Copyright 2018, Chris Coey and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using JuGP, JuMP
using Base.Test
import Pajarito, Ipopt, ECOS, GLPKMathProgInterface

mip_solver = GLPKMathProgInterface.GLPKSolverMIP(msg_lev=GLPK.MSG_OFF, tol_int=1e-9, tol_bnd=1e-8, mip_gap=1e-9)
nlp_solver = Ipopt.IpoptSolver(print_level=0)
conic_solver = ECOS.ECOSSolver(verbose=false)
intnlp_solver = Pajarito.PajaritoSolver(mip_solver=mip_solver, cont_solver=nlp_solver, log_level=0, rel_gap=1e-6)
intconic_solver = Pajarito.PajaritoSolver(mip_solver=mip_solver, cont_solver=conic_solver, log_level=0, rel_gap=1e-6)

methods = [:LogSumExp, :Conic]
cont_solvers = Dict()
int_solvers = Dict()
cont_solvers[:LogSumExp] = [nlp_solver]
cont_solvers[:Conic] = [conic_solver]
int_solvers[:LogSumExp] = [intnlp_solver]
int_solvers[:Conic] = [intconic_solver]

include("operators.jl")
include("models.jl")
include("applications.jl")
