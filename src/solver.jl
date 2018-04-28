#  Copyright 2016, Miles Lubin and contributors
#  Copyright 2018, Chris Coey and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

function solvehook(m::Model; suppress_warnings=false)
    gp = m.ext[:GP]::GPData
    m.solver.discretevalues = gp.discretevalues # hack
    solve(m, suppress_warnings=suppress_warnings, ignore_solve_hook=true)
end

mutable struct GPSolver <: MathProgBase.AbstractMathProgSolver
    method::Symbol
    real_solver
    discretevalues::Dict{Int,Vector{Float64}} # hack to pass this to "solver"
end

GPSolver(method, real_solver) = GPSolver(method, real_solver, Dict{Int,Vector{Float64}}())

mutable struct GPInternalModel <: MathProgBase.AbstractNonlinearModel
    status::Symbol
    method::Symbol # :LogSumExp or :Conic
    hasobj::Bool
    objinverted::Bool
    objconstant::Float64
    jump_model::JuMP.Model
    convexjl_problem::Convex.Problem
    real_solver
    y::Vector{JuMP.Variable}
    convexjl_y::Convex.Variable
    discretevalues::Dict{Int,Vector{Float64}}

    function GPInternalModel(method, real_solver, discretevalues)
        m = new(:Unsolved, method)
        m.real_solver = real_solver
        m.discretevalues = discretevalues
        return m
    end
end

MathProgBase.NonlinearModel(s::GPSolver) = GPInternalModel(s.method, s.real_solver, s.discretevalues)

using_jump(m::GPInternalModel) = (m.method == :LogSumExp)

function MathProgBase.loadproblem!(m::GPInternalModel, numvar, numconstr, x_lb, x_ub, g_lb, g_ub, sense, d::MathProgBase.AbstractNLPEvaluator)
    MathProgBase.initialize(d, [:ExprGraph])

    obj = nothing
    m.objinverted = false
    objorig = MathProgBase.obj_expr(d)
    if isa(objorig, Real)
        m.objconstant = objorig
        m.hasobj = false
    else
        objcheck = check_expr_gp(objorig)
        if isa(objcheck, Real)
            m.objconstant = objcheck
            m.hasobj = false
        else
            (obj, m.objconstant) = extract_constants(objcheck)
            if sense == :Max
                if !isa(obj, Monomial)
                    error("Can only maximize a monomial objective function")
                end
                obj = 1/obj
                m.objinverted = true
            end
            m.hasobj = true
        end
    end

    cons = Any[]
    cons_rhs = Float64[]

    for c in 1:numconstr
        con_expr = MathProgBase.constr_expr(d, c)
        # println("Constraint expr: ", con_expr)
        # Get constraint type
        @assert isexpr(con_expr, :call)
        @assert length(con_expr.args) == 3

        con_type = con_expr.args[1]
        # println("Constraint type: ", con_type)
        if con_type == :(>=) # flip
            con = check_expr_gp(Expr(:call, :*, -1, con_expr.args[2]))
        else
            con = check_expr_gp(con_expr.args[2])
        end

        if isa(con, Posynomial)
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
                    push!(new_con, mon)
                end
            end

            if found_negative
                # pull negative term to the right-hand side, then divide
                con = Posynomial(new_con)*(-1/negative_monomial) - 1
                #println(con)
            end
        end

        (con, con_constant) = extract_constants(con)
        push!(cons, con)
        if con_type == :(>=)
            push!(cons_rhs, (-g_lb[c] - con_constant))
        else
            push!(cons_rhs, (g_ub[c] - con_constant))
        end
    end

    if m.method == :LogSumExp
        if m.real_solver !== nothing
            jump_model = JuMP.Model(solver=m.real_solver)
        else
            jump_model = JuMP.Model()
        end

        # x = exp(y), y = log(x)
        y = @variable(jump_model, [1:numvar])
        for i in 1:numvar
            if x_lb[i] > 0
                setlowerbound(y[i], log(x_lb[i]))
            end
            @assert x_ub[i] > 0
            if isfinite(x_ub[i])
                setupperbound(y[i], log(x_ub[i]))
            end
        end

        if m.hasobj
            @objective(jump_model, Min, generate_epigraph(jump_model, y, obj))
        end

        for c in 1:numconstr
            con = cons[c]
            con_rhs = cons_rhs[c]
            @assert con_rhs > 0
            if g_lb[c] == g_ub[c] # equality constraint
                if isa(con, Posynomial)
                    @assert length(con.mons) == 1
                    con = con.mons[1]
                end
                @assert isa(con, Monomial)
                epi = generate_epigraph(jump_model, y, con)
                setlowerbound(epi, log(con_rhs))
                setupperbound(epi, log(con_rhs))
            else # <= or flipped >=
                epi = generate_epigraph(jump_model, y, con)
                setupperbound(epi, log(con_rhs))
            end
        end

        m.jump_model = jump_model
        m.y = y
    else
        # generate a conic model using Convex.jl
        @assert m.method == :Conic
        y_convex = Convex.Variable(numvar)

        if !m.hasobj
            prob = Convex.satisfy()
        else
            (extra_constraints, expr) = generate_epigraph(y_convex, obj)
            prob = Convex.minimize(expr, extra_constraints...)
        end

        for i in 1:numvar
            if x_lb[i] > 0
                push!(prob.constraints, y_convex[i] >= log(x_lb[i]))
            end
            @assert x_ub[i] > 0
            if isfinite(x_ub[i])
                push!(prob.constraints, y_convex[i] <= log(x_ub[i]))
            end
        end

        for c in 1:numconstr
            con = cons[c]
            con_rhs = cons_rhs[c]
            @assert con_rhs > 0
            if g_lb[c] == g_ub[c] # equality constraint
                if isa(con, Posynomial)
                    @assert length(con.mons) == 1
                    con = con.mons[1]
                end
                @assert isa(con, Monomial)
                extra_constraints, expr = generate_epigraph(y_convex, con)
                @assert length(extra_constraints) == 0
                push!(prob.constraints, expr == log(con_rhs))
            else # <= or flipped >=
                extra_constraints, expr = generate_epigraph(y_convex, con)
                append!(prob.constraints, extra_constraints)
                push!(prob.constraints, expr <= log(con_rhs))
            end
        end

        m.convexjl_problem = prob
        m.convexjl_y = y_convex
    end

    # generate the formulation for discrete part of the model
    discretevalues!(m)
end

function discretevalues!(m::GPInternalModel)
    for (ind, values) in m.discretevalues
        N = length(values)
        @assert all(v -> (v > 0), values)
        logvals = log(values)

        if using_jump(m)
            # dummy formulation
            discrete = @variable(m.jump_model, [1:N], category=:Bin)
            @constraint(m.jump_model, sum(discrete) == 1)
            @constraint(m.jump_model, m.y[ind] == dot(logvals, discrete))
        else
            discrete = Convex.Variable(N, :Bin)
            push!(m.convexjl_problem.constraints, sum(discrete) == 1)
            push!(m.convexjl_problem.constraints, m.convexjl_y[ind] == dot(logvals, discrete))
        end
    end
end

function generate_aux(m::JuMP.Model, y, mon::Monomial)
    # generate an auxiliary variable representing the linear part of the monomial in log space
    @assert mon.c > 0
    aux = @variable(m)
    @constraint(m, aux == log(mon.c) + sum(v*y[i] for (i,v) in mon.terms))
    return aux
end

function generate_epigraph(m::JuMP.Model, y, mon::Monomial)
    return generate_aux(m , y, mon)
end

function generate_epigraph(m::JuMP.Model, y, pos::Posynomial)
    # generate an auxiliary variable representing the logsumexp of the given posynomial
    if length(pos.mons) == 1
        return generate_aux(m, y, pos.mons[1])
    else
        epi = @variable(m)
        mon_aux = [generate_aux(m, y, mon) for mon in pos.mons]
        # could also generate the separable version here
        @NLconstraint(m, log(sum(exp(aux) for aux in mon_aux)) <= epi)
        return epi
    end
end

# returns a set of constraints to add as well
function generate_epigraph(y::Convex.Variable, mon::Monomial)
    return ([], log(mon.c) + sum(v*y[i] for (i,v) in mon.terms))
end

function generate_epigraph(y::Convex.Variable, pos::Posynomial)
    if length(pos.mons) == 1
        return generate_epigraph(y, pos.mons[1])
    else
        aux = Convex.Variable(length(pos.mons))
        cons = [aux[k] == generate_epigraph(y, pos.mons[k])[2] for k in 1:length(pos.mons)]
        return (cons, Convex.logsumexp(aux))
    end
end

MathProgBase.setwarmstart!(m::GPInternalModel, x) = nothing

function MathProgBase.optimize!(m::GPInternalModel)
    if using_jump(m)
        m.status = solve(m.jump_model)
    else
        if m.real_solver !== nothing
            Convex.solve!(m.convexjl_problem, m.real_solver)
        else
            Convex.solve!(m.convexjl_problem)
        end
        m.status = m.convexjl_problem.status
    end
    return nothing
end

MathProgBase.status(m::GPInternalModel) = m.status

function MathProgBase.getobjval(m::GPInternalModel)
    if m.hasobj
        if using_jump(m)
            val = exp(getobjectivevalue(m.jump_model))
        else
            val = exp(m.convexjl_problem.optval)
        end
        if m.objinverted
            val = 1/val
        end
    else
        val = 0.
    end
    val += m.objconstant
    return val
end

function MathProgBase.getsolution(m::GPInternalModel)
    if using_jump(m)
        y = getvalue(m.y)
    else
        y = Convex.evaluate(m.convexjl_y)
        if isa(y, Float64) # bizarre behavior
            y = [y]
        else
            y = vec(y)
        end
    end
    return exp.(y)
end
