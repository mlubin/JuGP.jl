#  Copyright 2016, Miles Lubin and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

function solvehook(m::Model; suppress_warnings=false)

    gp = m.ext[:GP]::GPData
    # hack!
    m.solver.discretevalues = gp.discretevalues

    solve(m,suppress_warnings=suppress_warnings,ignore_solve_hook=true)
end



type GPSolver <: MathProgBase.AbstractMathProgSolver
    real_solver
    discretevalues::Dict{Int,Vector{Float64}} # hack to pass this to "solver"
end
GPSolver(real_solver) = GPSolver(real_solver,Dict{Int,Vector{Float64}}())

type GPInternalModel <: MathProgBase.AbstractNonlinearModel
    status::Symbol
    jump_model::JuMP.Model
    real_solver
    y::Vector{JuMP.Variable}
    discretevalues::Dict{Int,Vector{Float64}}
    function GPInternalModel(real_solver,discretevalues)
        m = new(:Unsolved)
        m.real_solver = real_solver
        m.discretevalues = discretevalues
        return m
    end
end
MathProgBase.NonlinearModel(s::GPSolver) = GPInternalModel(s.real_solver,s.discretevalues)

function MathProgBase.loadproblem!(m::GPInternalModel, numVar, numConstr, x_lb, x_ub, g_lb, g_ub, sense,
    d::MathProgBase.AbstractNLPEvaluator)

    MathProgBase.initialize(d, [:ExprGraph])

    obj = check_expr_gp(MathProgBase.obj_expr(d))
    if sense == :Max
        @assert isa(obj, Monomial)
        obj = 1/obj
    end
    obj, obj_constant = extract_constants(obj)
    @assert obj_constant == 0.0

    cons = Any[]
    cons_rhs = Float64[]
    
    for c in 1:numConstr
        con_expr = MathProgBase.constr_expr(d,c)
        # Get constraint type
        @assert isexpr(con_expr,:comparison)
        @assert length(con_expr.args) == 3
  
        con_type = con_expr.args[2]
        #println("Constraint type: $con_type")
        if con_type == :(>=) # flip
            con = check_expr_gp(Expr(:call,:*,-1,con_expr.args[1]))
        else
            con = check_expr_gp(con_expr.args[1])
        end

        if isa(con,Posynomial)
            # check if any monomials have a negative coefficient, we can have at most one
            found_negative = false
            negative_monomial = Monomial(1.0)
            new_con = []
            for mon in con.mons
                if mon.c < 0
                    @assert !found_negative
                    found_negative = true
                    negative_monomial = mon
                else
                    push!(new_con,mon)
                end
            end
            if found_negative
                # pull negative term to the right-hand side, then divide
                con = Posynomial(new_con)*(-1/negative_monomial)-1
                #println(con)
            end
        end
        con, con_constant = extract_constants(con)
        push!(cons, con)
        if con_type == :(>=)
            push!(cons_rhs, -g_lb[c]-con_constant)
        else
            push!(cons_rhs, g_ub[c]-con_constant)
        end

    end

    if m.real_solver !== nothing
        jump_model = JuMP.Model(solver=m.real_solver)
    else
        jump_model = JuMP.Model()
    end
    # x = exp(y), y = log(x)
    @defVar(jump_model, y[1:numVar])
    for i in 1:numVar
        if x_lb[i] > 0
            setLower(y[i], log(x_lb[i]))
        end
        @assert x_ub[i] > 0
        if isfinite(x_ub[i])
            setUpper(y[i], log(x_ub[i]))
        end
    end

    @setObjective(jump_model, Min, generate_epigraph(jump_model, y, obj))

    for c in 1:numConstr
        con = cons[c]
        con_rhs = cons_rhs[c]
        #@show con_rhs
        @assert con_rhs > 0
        if g_lb[c] == g_ub[c] # equality constraint
            if isa(con, Posynomial)
                @assert length(con.mons) == 1
                con = con.mons[1]
            end
            @assert isa(con, Monomial)
            epi = generate_epigraph(jump_model, y, con)
            setLower(epi, log(con_rhs))
            setUpper(epi, log(con_rhs))
        else # <= or flipped >=
            epi = generate_epigraph(jump_model, y, con)
            setUpper(epi, log(con_rhs))
        end
    end

    m.jump_model = jump_model
    m.y = y

    # generate the formulation for discrete part of the model
    discretevalues!(m)

end

function discretevalues!(m::GPInternalModel)
    for (ind, values) in m.discretevalues
        N = length(values)
        @assert all(v -> v > 0, values)
        logvals = log(values)
        # dummy formulation
        @defVar(m.jump_model, discrete[1:N], Bin)
        @addConstraint(m.jump_model, sum(discrete) == 1)
        @addConstraint(m.jump_model, m.y[ind] == dot(logvals,discrete))
    end
end

function generate_aux(m::JuMP.Model, y, mon::Monomial)
    # generate an auxiliary variable representing
    # the linear part of the monomial in log space.
    @assert mon.c > 0

    @defVar(m, aux)
    @addConstraint(m, aux == log(mon.c) + sum{v*y[i], (i,v) in mon.terms})
    return aux
end

function generate_epigraph(m::JuMP.Model, y, mon::Monomial)
    return generate_aux(m , y, mon)
end

function generate_epigraph(m::JuMP.Model, y, pos::Posynomial)
    # generate an auxiliary variable representing the logsumexp
    # of the given posynomial

    if length(pos.mons) == 1
        return generate_aux(m, y, pos.mons[1])
    else
        @defVar(m, epi)
        mon_aux = [generate_aux(m, y, mon) for mon in pos.mons]
        # could also generate the separable version here
        @addNLConstraint(m, log(sum{exp(aux), aux in mon_aux}) <= epi)
        return epi
    end
end

MathProgBase.setwarmstart!(m::GPInternalModel,x) = nothing
function MathProgBase.optimize!(m::GPInternalModel)
    m.status = solve(m.jump_model)
end
MathProgBase.status(m::GPInternalModel) = m.status

function MathProgBase.getobjval(m::GPInternalModel)
    return exp(getObjectiveValue(m.jump_model))
end

function MathProgBase.getsolution(m::GPInternalModel)
    y = getValue(m.y)
    return exp(y)
end
