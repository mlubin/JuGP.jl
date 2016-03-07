using JuGP, JuMP

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

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

@testset "Optimize the shape of a box" begin
    #=
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

    # comparing optimal answers with http://gpkit.readthedocs.org/en/latest/examples.html#maximizing-the-volume-of-a-box
    @test_approx_eq_eps getValue(d) 8.17 1e-2
    @test_approx_eq_eps getValue(h) 8.163 1e-2
    @test_approx_eq_eps getValue(w) 4.081 1e-2
    # this actually should be inverted
    @test_approx_eq_eps getObjectiveValue(m) 0.003674 1e-2
end
