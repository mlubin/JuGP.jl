using JuGP, JuMP

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

include("operators.jl")

@testset "Model tests" begin
@testset "Equality constraints" begin
    m = Model(solver=GPSolver())

    @defVar(m, x)
    @setNLObjective(m, Min, x)
    @addNLConstraint(m, 2x == 4)

    status = solve(m)
    @test status == :Optimal
    @test_approx_eq_eps getValue(x) 2.0 1e-4
end

@testset "Monomial^Number" begin
    m = Model(solver=GPSolver())

    @defVar(m, x)
    @setNLObjective(m, Min, x)
    @addNLConstraint(m, x^2 == 4)

    status = solve(m)
    @test status == :Optimal
    @test_approx_eq_eps getValue(x) 2.0 1e-4
end

@testset "Monomial+Number and Number+Monomial" begin
    m = Model(solver=GPSolver())

    @defVar(m, x)
    @setNLObjective(m, Min, x)
    @addNLConstraint(m, x^2 + 6 == 10)

    status = solve(m)
    @test status == :Optimal
    @test_approx_eq_eps getValue(x) 2.0 1e-4

    m = Model(solver=GPSolver())

    @defVar(m, x)
    @setNLObjective(m, Min, x)
    @addNLConstraint(m, 6 + x^2 == 10)

    status = solve(m)
    @test status == :Optimal
    @test_approx_eq_eps getValue(x) 2.0 1e-4
end

@testset "Optimize the shape of a box" begin
    #=
    From http://gpkit.readthedocs.org/en/latest/examples.html#maximizing-the-volume-of-a-box

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

    m = Model(solver=GPSolver())

    @defVar(m, h)
    @defVar(m, w)
    @defVar(m, d)

    @setNLObjective(m, Max, h*w*d)

    @addNLConstraint(m, 2(h*w+h*d) ≤ Awall)
    @addNLConstraint(m, w*d ≤ Afloor)

    @addNLConstraint(m, α ≤ h/w)
    @addNLConstraint(m, h/w ≤ β)

    @addNLConstraint(m, γ ≤ d/w)
    @addNLConstraint(m, d/w ≤ δ)

    status = solve(m)

    @test status == :Optimal

    # comparing optimal answers with source
    @test_approx_eq_eps getValue(d) 8.17 1e-2
    @test_approx_eq_eps getValue(h) 8.163 1e-2
    @test_approx_eq_eps getValue(w) 4.081 1e-2
    # this actually should be inverted
    @test_approx_eq_eps getObjectiveValue(m) 0.003674 1e-2

    # with >= constraints
    m = Model(solver=GPSolver())

    @defVar(m, h)
    @defVar(m, w)
    @defVar(m, d)

    @setNLObjective(m, Max, h*w*d)

    @addNLConstraint(m, Awall ≥ 2(h*w+h*d))
    @addNLConstraint(m, Afloor ≥ w*d)

    @addNLConstraint(m, h/w ≥ α)
    @addNLConstraint(m, β ≥ h/w)

    @addNLConstraint(m, d/w ≥ γ)
    @addNLConstraint(m, δ ≥ d/w)

    status = solve(m)

    @test status == :Optimal

    @test_approx_eq_eps getValue(d) 8.17 1e-2
    @test_approx_eq_eps getValue(h) 8.163 1e-2
    @test_approx_eq_eps getValue(w) 4.081 1e-2
    # this actually should be inverted
    @test_approx_eq_eps getObjectiveValue(m) 0.003674 1e-2
end

@testset "Optimize gate sizes (optional integer variables)" begin
    #=
    See Boyd GP tutorial 2007, "gate sizing" examples
    Translated fom YALMIP formulation from
    http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Tutorials.GeometricProgramming
    =#

    a     = ones(7)
    alpha = ones(7)
    beta  = ones(7)
    gamma = ones(7)
    f = [1, 0.8, 1, 0.7, 0.7, 0.5, 0.5]
    e = [1, 2, 1, 1.5, 1.5, 1, 2]
    Cout6 = 10
    Cout7 = 10

    m = Model(solver=GPSolver())

    @defVar(m, D)

    @defVar(m, x[i=1:7] >= 1)

    @setNLObjective(m, Min, D)

    @defNLExpr(m, C[i=1:7], alpha[i] + beta[i] * x[i])
    @defNLExpr(m, A, sum{a[i] * x[i], i=1:7})
    @defNLExpr(m, P, sum{f[i] * e[i] * x[i], i=1:7})
    @defNLExpr(m, R[i=1:7], gamma[i] / x[i])
    @defNLExpr(m, D1, R[1] * C[4])
    @defNLExpr(m, D2, R[2] * (C[4] + C[5]))
    @defNLExpr(m, D3, R[3] * (C[5] + C[7]))
    @defNLExpr(m, D4, R[4] * (C[6] + C[7]))
    @defNLExpr(m, D5, R[5] * C[7])
    @defNLExpr(m, D6, R[6] * Cout6)
    @defNLExpr(m, D7, R[7] * Cout7)

    @addNLConstraint(m, P <= 20)
    @addNLConstraint(m, A <= 100)
    @addNLConstraint(m, D1+D4+D6 <= D)
    @addNLConstraint(m, D1+D4+D7 <= D)
    @addNLConstraint(m, D2+D4+D6 <= D)
    @addNLConstraint(m, D2+D4+D7 <= D)
    @addNLConstraint(m, D2+D5+D7 <= D)
    @addNLConstraint(m, D3+D5+D6 <= D)
    @addNLConstraint(m, D3+D7 <= D)

    # test continuous problem
    status = solve(m)
    @test status == :Optimal
    # comparing optimal continuous answers with YALMIP solutions
    x_opt = [1.9563, 3.1588, 3.0455, 3.3454, 1.6713, 3.1224, 3.1155]
    @test_approx_eq_eps getValue(x) x_opt 1e-3
    @test_approx_eq_eps getObjectiveValue(m) 7.8936 1e-4

    # test integer-constrained problem
    # setCategory(x, Int)
    # status = solve(m)
    # @test status == :Optimal
    # # comparing optimal continuous answers with YALMIP solutions
    # x_opt = [1.9563, 3.1588, 3.0455, 3.3454, 1.6713, 3.1224, 3.1155]
    # @test_approx_eq_eps getValue(x) x_opt 1e-4
    # @test_approx_eq_eps getObjectiveValue(m) 7.8936 1e-4
end

include("cvx_examples.jl")
end
