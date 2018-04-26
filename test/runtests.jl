using JuGP, JuMP
using Base.Test

include(Pkg.dir("JuMP", "test", "solvers.jl"))

methods = [:LogSumExp, :Conic]
cont_solvers = Dict()
cont_solvers[:LogSumExp] = []
cont_solvers[:Conic] = []

# ipt && push!(cont_solvers[:LogSumExp], Ipopt.IpoptSolver(print_level=0))
eco && push!(cont_solvers[:Conic], ECOS.ECOSSolver(verbose=false))

include("operators.jl")

@testset "Model tests" begin
    @testset "Equality constraints, $meth" for meth in methods
        for solv in cont_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, x)
            @NLobjective(m, Min, x)
            @NLconstraint(m, 2x == 4)

            status = solve(m)
            @test status == :Optimal
            @test getvalue(x) ≈ 2.0 atol=1e-4
            @test getobjectivevalue(m) ≈ 0.804067384844554 atol=1.0e-6
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
            @test getvalue(x) ≈ 2.0 atol=1e-4
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
            @test getvalue(x) ≈ 2.0 atol=1e-4

            m = GPModel(method=meth, solver=solv)

            @variable(m, x)
            @NLobjective(m, Min, x)
            @NLconstraint(m, 6 + x^2 == 10)

            status = solve(m)
            @test status == :Optimal
            @test getvalue(x) ≈ 2.0 atol=1e-4
        end
    end

    @testset "Optimize the shape of a box, $meth" for meth in methods
        #=
        From http://gpkit.readthedocs.org/en/latest/examples.html#maximizing-the-volume-of-a-box
        We copy the license for gpkit below:
        "Copyright (c) 2015 MIT Convex Optimization Group

        Permission is hereby granted, free of charge, to any person obtaining a copy
        of this software and associated documentation files (the "Software"), to deal
        in the Software without restriction, including without limitation the rights
        to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
        copies of the Software, and to permit persons to whom the Software is
        furnished to do so, subject to the following conditions:

        The above copyright notice and this permission notice shall be included in all
        copies or substantial portions of the Software.""

        Optimize shape of a box-shaped structure with
        - height h
        - width w
        - depth d
        Limit on wall area: 2(hw+hd)
        Limit on floor area: wd
        Bounds on ratios: h/w, w/d
        max hwd
         st 2(hw + hd) ≤ Awall
            wd ≤ Afloor
            α ≤ h/w ≤ β
            γ ≤ d/w ≤ δ

        or in standard form
        min h⁻¹ w⁻¹ d⁻¹
         st (2/Awall)*hw + (2/Awall)*hd ≤ 1
            (1/Aflr)wd ≤ 1
            αh⁻¹w ≤ 1
            (1/β)hw⁻¹ ≤ 1
            γw⁻¹d ≤ 1
            (1/δ)w⁻¹d ≤ 1
        =#

        Afloor = 50
        Awall  = 200
        α = 2
        β = 10
        γ = 2
        δ = 10

        for solv in cont_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, h)
            @variable(m, w)
            @variable(m, d)

            @NLobjective(m, Max, h*w*d)

            @NLconstraint(m, 2(h*w+h*d) ≤ Awall)
            @NLconstraint(m, w*d ≤ Afloor)

            @NLconstraint(m, α ≤ h/w)
            @NLconstraint(m, h/w ≤ β)

            @NLconstraint(m, γ ≤ d/w)
            @NLconstraint(m, d/w ≤ δ)

            status = solve(m)

            @test status == :Optimal

            # comparing optimal answers with source
            @test getvalue(d) ≈ 8.17 atol=1e-2
            @test getvalue(h) ≈ 8.163 atol=1e-2
            @test getvalue(w) ≈ 4.081 atol=1e-2
            # this actually should be inverted
            @test getobjectivevalue(m) ≈ 0.003674 atol=1e-2

            # with >= constraints
            m = GPModel(method=meth, solver=solv)

            @variable(m, h)
            @variable(m, w)
            @variable(m, d)

            @NLobjective(m, Max, h*w*d)

            @NLconstraint(m, Awall ≥ 2(h*w+h*d))
            @NLconstraint(m, Afloor ≥ w*d)

            @NLconstraint(m, h/w ≥ α)
            @NLconstraint(m, β ≥ h/w)

            @NLconstraint(m, d/w ≥ γ)
            @NLconstraint(m, δ ≥ d/w)

            status = solve(m)

            @test status == :Optimal

            @test getvalue(d) ≈ 8.17 atol=1e-2
            @test getvalue(h) ≈ 8.163 atol=1e-2
            @test getvalue(w) ≈ 4.081 atol=1e-2
            # this actually should be inverted
            @test getobjectivevalue(m) ≈ 0.003674 atol=1e-2
        end
    end

    if nlw && osl # these are defined in JuMP's test/solvers.jl
        # Only run if Bonmin is installed
        @testset "Optimize gate sizes (optional integer variables)" begin
            #=
            See Boyd GP tutorial 2007, "gate sizing" examples
            Translated from YALMIP formulation:
            http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Tutorials.GeometricProgramming
            Received permission from Johan Löfberg to distribute under MPL
            =#

            a     = ones(7)
            alpha = ones(7)
            beta  = ones(7)
            gamma = ones(7)
            f = [1, 0.8, 1, 0.7, 0.7, 0.5, 0.5]
            e = [1, 2, 1, 1.5, 1.5, 1, 2]
            Cout6 = 10
            Cout7 = 10

            m = GPModel(solver=AmplNLWriter.BonminNLSolver())

            @variable(m, D)

            @variable(m, x[i=1:7] >= 1)

            @NLobjective(m, Min, D)

            @NLexpression(m, C[i=1:7], alpha[i] + beta[i] * x[i])
            @NLexpression(m, A, sum(a[i] * x[i] for i=1:7))
            @NLexpression(m, P, sum(f[i] * e[i] * x[i] for i=1:7))
            @NLexpression(m, R[i=1:7], gamma[i] / x[i])
            @NLexpression(m, D1, R[1] * C[4])
            @NLexpression(m, D2, R[2] * (C[4] + C[5]))
            @NLexpression(m, D3, R[3] * (C[5] + C[7]))
            @NLexpression(m, D4, R[4] * (C[6] + C[7]))
            @NLexpression(m, D5, R[5] * C[7])
            @NLexpression(m, D6, R[6] * Cout6)
            @NLexpression(m, D7, R[7] * Cout7)

            @NLconstraint(m, P <= 20)
            @NLconstraint(m, A <= 100)
            @NLconstraint(m, D1+D4+D6 <= D)
            @NLconstraint(m, D1+D4+D7 <= D)
            @NLconstraint(m, D2+D4+D6 <= D)
            @NLconstraint(m, D2+D4+D7 <= D)
            @NLconstraint(m, D2+D5+D7 <= D)
            @NLconstraint(m, D3+D5+D6 <= D)
            @NLconstraint(m, D3+D7 <= D)

            # test continuous problem
            status = solve(m)
            @test status == :Optimal
            # comparing optimal continuous answers with YALMIP solutions
            x_opt = [1.9563, 3.1588, 3.0455, 3.3454, 1.6713, 3.1224, 3.1155]
            @test getvalue(x) ≈ x_opt atol=1e-3
            @test getobjectivevalue(m) ≈ 7.8936 atol=1e-4

            # test integer-constrained problem
            for i in 1:7
                setDiscrete(x[i], 1:10)
            end
            status = solve(m)
            @test status == :Optimal
            # comparing optimal integer answers with YALMIP solutions
            x_opt = [2, 3, 3, 3, 2, 3, 3]
            @test getvalue(x) ≈ x_opt atol=1e-4
            @test getobjectivevalue(m) ≈ 8.3333 atol=1e-4
        end
    end
end

include("cvx_examples.jl")
