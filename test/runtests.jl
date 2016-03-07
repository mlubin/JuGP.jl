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

@testset "Optimal doping profile" begin
    # Problem from Boyd's GP tutorial
    # Based on public domain CVX example:
    # http://web.cvxr.com/cvx/examples/gp_tutorial/html/basic_odp.html

    # Discretization size
    M = 50
    # Problem constants
    g1 = 0.42
    g2 = 0.69
    Nmax = 5*10.0^18
    Nmin = 5*10.0^16
    Nref = 10.0^17
    Dn0  = 20.72
    ni0  = 1.4*10.0^10
    WB   = 10.0^-5
    C    =  WB^2/((M^2)*(Nref^g1)*Dn0)
    # Exponent powers
    pwi  = g2 -1;
    pwj  = 1+g1-g2;

    doping = Model(solver=GPSolver())
    @defVar(doping, v[1:M])
    @defVar(doping, y[1:M])
    @defVar(doping, w[1:M])

    # Minimize base transmit times
    @setNLObjective(doping, Min, C * w[1])

    for i in 1:M-1
        @addNLConstraint(doping, y[i+1] + v[i]^pwj <= y[i])
        @addNLConstraint(doping, w[i+1] + y[i]*v[i]^pwi <= w[i])
    end
    @addNLConstraint(doping, y[M] == v[M]^pwj)
    @addNLConstraint(doping, w[M] == y[M]*v[M]^pwi)

    solve(doping)
end
