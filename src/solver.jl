#  Copyright 2016, Miles Lubin and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

export GPSolver
type GPSolver <: MathProgBase.AbstractMathProgSolver
end
type GPModel <: MathProgBase.AbstractNonlinearModel
    status::Symbol
    jump_model::JuMP.Model
    y::Vector{JuMP.Variable}
    function GPModel()
        return new(:Unsolved)
    end
end
MathProgBase.NonlinearModel(s::GPSolver) = GPModel()

function MathProgBase.loadproblem!(m::GPModel, numVar, numConstr, x_lb, x_ub, g_lb, g_ub, sense, 
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

    jump_model = JuMP.Model()
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

MathProgBase.setwarmstart!(m::GPModel,x) = nothing
function MathProgBase.optimize!(m::GPModel)
    m.status = solve(m.jump_model)
end
MathProgBase.status(m::GPModel) = m.status

function MathProgBase.getobjval(m::GPModel)
    return exp(getObjectiveValue(m.jump_model))
end

function MathProgBase.getsolution(m::GPModel)
    y = getValue(m.y)
    return exp(y)
end
