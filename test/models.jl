#  Copyright 2016, Miles Lubin and contributors
#  Copyright 2018, Chris Coey and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using JuGP, JuMP
import JuGP: Monomial, Posynomial
using Base.Test

@testset "Model tests" begin
    @testset "LP model, $meth" for meth in methods
        for solv in cont_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, 2 <= x <= 2)
            @variable(m, y >= 1.5)
            @objective(m, Min, y)
            @constraint(m, x == y)

            status = solve(m)
            @test status == :Optimal
            @test getvalue(x) ≈ 2 atol=1e-6
            @test getvalue(x) ≈ 2 atol=1e-6
            @test getobjectivevalue(m) ≈ 2 atol=1e-6
        end
    end

    @testset "QP objective model, $meth" for meth in methods
        for solv in cont_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, 2 <= x <= 2)
            @variable(m, y >= 1.5)
            @objective(m, Min, x*y)
            @constraint(m, x <= y)

            status = solve(m)
            @test status == :Optimal
            @test getvalue(x) ≈ 2 atol=1e-6
            @test getvalue(x) ≈ 2 atol=1e-6
            @test getobjectivevalue(m) ≈ 4 atol=1e-6
        end
    end

    @testset "QP constraint model, $meth" for meth in methods
        for solv in cont_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, 2 <= x <= 2)
            @variable(m, y >= 1.5)
            @objective(m, Min, y)
            @constraint(m, x*y == 4)

            status = solve(m)
            @test status == :Optimal
            @test getvalue(x) ≈ 2 atol=1e-6
            @test getvalue(x) ≈ 2 atol=1e-6
            @test getobjectivevalue(m) ≈ 2 atol=1e-6
        end
    end

    @testset "Equality constraints, $meth" for meth in methods
        for solv in cont_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, x)
            @NLobjective(m, Min, x)
            @NLconstraint(m, 2x == 4)
            @constraint(m, x == 2)

            status = solve(m)
            @test status == :Optimal
            @test getvalue(x) ≈ 2 atol=1e-6
            @test getobjectivevalue(m) ≈ 2 atol=1e-6
        end
    end

    @testset "Inequality constraints, $meth" for meth in methods
        for solv in cont_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, x)
            @objective(m, Min, x)
            @NLconstraint(m, 2x >= 4)
            @constraint(m, x <= 2)

            status = solve(m)
            @test status == :Optimal
            @test getvalue(x) ≈ 2 atol=1e-6
            @test getobjectivevalue(m) ≈ 2 atol=1e-6
        end
    end

    @testset "Monomial^Number, $meth" for meth in methods
        for solv in cont_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, x)
            @NLobjective(m, Min, x)
            @NLconstraint(m, x^2 == 4)

            status = solve(m)
            @test status == :Optimal
            @test getvalue(x) ≈ 2 atol=1e-4
            @test getobjectivevalue(m) ≈ 2 atol=1e-4
        end
    end

    @testset "Monomial+Number and Number+Monomial, $meth" for meth in methods
        for solv in cont_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, x)
            @NLobjective(m, Min, x)
            @NLconstraint(m, x^2 + 6 == 10)

            status = solve(m)
            @test status == :Optimal
            @test getvalue(x) ≈ 2 atol=1e-4
            @test getobjectivevalue(m) ≈ 2 atol=1e-4

            m = GPModel(method=meth, solver=solv)

            @variable(m, x)
            @NLobjective(m, Min, x)
            @NLconstraint(m, 6 + x^2 == 10)

            status = solve(m)
            @test status == :Optimal
            @test getvalue(x) ≈ 2 atol=1e-4
            @test getobjectivevalue(m) ≈ 2 atol=1e-4
        end
    end

    @testset "Discrete values, $meth" for meth in methods
        for solv in int_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, x)
            @objective(m, Min, x)
            @NLconstraint(m, x^2 >= 4.5)

            setdiscretevalues(x, 1:5)

            status = solve(m)
            @test status == :Optimal
            @test getvalue(x) ≈ 3 atol=1e-6
            @test getobjectivevalue(m) ≈ 3 atol=1e-6
        end
    end

    @testset "Monomial objective, $meth" for meth in methods
        for solv in cont_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, x >= 2)
            @variable(m, y >= 1.5)
            @NLobjective(m, Min, x*y)
            @constraint(m, x == y)

            status = solve(m)
            @test status == :Optimal
            @test getvalue(x) ≈ 2 atol=1e-6
            @test getvalue(y) ≈ 2 atol=1e-6
            @test getobjectivevalue(m) ≈ 4 atol=1e-6
        end
    end

    @testset "Posynomial objective, $meth" for meth in methods
        for solv in cont_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, x >= 2)
            @variable(m, y >= 1.5)
            @NLobjective(m, Min, x*y + (x + y)/2)
            @NLconstraint(m, x == y)

            status = solve(m)
            @test status == :Optimal
            @test getvalue(x) ≈ 2 atol=1e-6
            @test getvalue(y) ≈ 2 atol=1e-6
            @test getobjectivevalue(m) ≈ 6 atol=1e-6
        end
    end

    @testset "Max objective, $meth" for meth in methods
        for solv in cont_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, x >= 2)
            @variable(m, y >= 1.5)
            @NLobjective(m, Max, 1/(x*y))
            @constraint(m, x == y)

            status = solve(m)
            @test status == :Optimal
            @test getvalue(x) ≈ 2 atol=1e-6
            @test getvalue(y) ≈ 2 atol=1e-6
            @test getobjectivevalue(m) ≈ 1/4 atol=1e-6
        end
    end

    @testset "Objective with constant, $meth" for meth in methods
        for solv in cont_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, x >= 2)
            @variable(m, y >= 1.5)
            @NLobjective(m, Min, x*y - 2)
            @NLconstraint(m, x == y)

            status = solve(m)
            @test status == :Optimal
            @test getvalue(x) ≈ 2 atol=1e-6
            @test getvalue(y) ≈ 2 atol=1e-6
            @test getobjectivevalue(m) ≈ 2 atol=1e-6
        end
    end

    @testset "Max objective with constant, $meth" for meth in methods
        for solv in cont_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, x >= 2)
            @variable(m, y >= 1.5)
            @NLobjective(m, Max, 1/(x*y) + 3/4)
            @NLconstraint(m, x == y)

            status = solve(m)
            @test status == :Optimal
            @test getvalue(x) ≈ 2 atol=1e-6
            @test getvalue(x) ≈ 2 atol=1e-6
            @test getobjectivevalue(m) ≈ 1 atol=1e-6
        end
    end

    @testset "Constant objective, $meth" for meth in methods
        for solv in cont_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, 2 <= x <= 2)
            @variable(m, y >= 1.5)
            @NLobjective(m, Min, 1)
            @NLconstraint(m, x == y)

            status = solve(m)
            @test status == :Optimal
            @test getvalue(x) ≈ 2 atol=1e-6
            @test getvalue(x) ≈ 2 atol=1e-6
            @test getobjectivevalue(m) ≈ 1 atol=1e-6
        end
    end
end
