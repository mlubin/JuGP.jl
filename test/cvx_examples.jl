#  Copyright 2016, Miles Lubin and contributors
#  Copyright 2018, Chris Coey and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# These test cases are translated from the CVX example library
# cvxr.com/cvx/examples/
#
# Below, we copy the licence from cvxr.com/cvx/licensing/
# "The contents of the example library, which is distributed with CVX in the
# examples/ subdirectory, is public domain. You are free to use them in any
# way you wish; but when you do, we request that you give appropriate credit
# to the authors. A number of people have contributed to the examples in this
# library, including Lieven Vandenberghe, Joëlle Skaf, Argyris Zymnis, Almir
# Mutapcic, Michael Grant, and Stephen Boyd."

@testset "CVX examples" begin
    @testset "Optimal doping profiles, $meth" for meth in methods
        #=
        Boyd, Kim, Vandenberghe, and Hassibi, "A tutorial on geometric programming"
        Joshi, Boyd, and Dutton, "Optimal doping profiles via geometric programming"
        web.cvxr.com/cvx/examples/gp_tutorial/html/basic_odp.html

        Determines the optimal doping profile that minimizes base transit
        time in a (homojunction) bipolar junction transistor.
        This problem can be posed as a GP:

           minimize   tau_B
               s.t.   Nmin <= v <= Nmax
                      y_(i+1) + v_i^const1 <= y_i
                      w_(i+1) + v_i^const2 <= w_i, etc...

        where variables are v_i, y_i, and w_i.
        =#

        # discretization size
        M = 50

        # problem constants
        g1 = 0.42
        g2 = 0.69
        Nmax = 5*10^18
        Nmin = 5*10^16
        Nref = 10^17
        Dn0 = 20.72
        ni0 = 1.4*(10^10)
        WB = 10.0^(-5)
        C =  WB^2/((M^2)*(Nref^g1)*Dn0)

        # exponent powers
        pwi = g2-1
        pwj = 1+g1-g2

        for solv in cont_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, Nmin <= v[i=1:M] <= Nmax)
            @variable(m, y[i=1:M])
            @variable(m, w[i=1:M])

            @NLobjective(m, Min, C*w[1]) # the base transmit time

            @NLconstraint(m, ly[i=1:M-1], y[i+1] + v[i]^pwj <= y[i])
            @NLconstraint(m, lw[i=1:M-1], w[i+1] + y[i]*v[i]^pwi <= w[i])
            @NLconstraint(m, ey, y[M] == v[M]^pwj)
            @NLconstraint(m, ew, w[M] == y[M]*v[M]^pwi)

            # test continuous problem
            status = solve(m)

            @test status == :Optimal

            # comparE optimal continuous answers with cvx solutions
            @test getobjectivevalue(m) ≈ 1.57873e-12 atol=1e-5
            v_opt = [
                4.99999998094197e18,4.99999851343142e18,2.54007808655285e18,1.50011889500592e18,
                9.78811490079119e17,6.83192189106757e17,5.00671118987733e17,3.80719724163692e17,
                2.9802359179552e17,2.38808625851304e17,1.95080971252896e17,1.61954847488488e17,
                1.36314038747507e17,1.16098679267784e17,9.9905244235857e16,8.67520833800946e16,
                7.59368040516422e16,6.69466735606206e16,5.94007147840608e16,5.30113996473861e16,
                5.00000029431718e16,5.00000009149364e16,5.00000005457447e16,5.00000003896229e16,
                5.00000003032521e16,5.00000002483656e16,5.0000000210376e16,5.00000001825071e16,
                5.00000001611816e16,5.00000001443311e16,5.00000001306762e16,5.00000001193836e16,
                5.00000001098865e16,5.00000001017856e16,5.00000000947924e16,5.00000000886924e16,
                5.00000000833228e16,5.00000000785583e16,5.00000000743003e16,5.00000000704698e16,
                5.00000000670038e16,5.00000000638508e16,5.0000000060967e16,5.00000000583167e16,
                5.00000000558689e16,5.00000000535966e16,5.00000000514724e16,5.00000000494647e16,
                5.00000000475203e16,5.00000000458168e16
                ]
            @test log.(getvalue(v)) ≈ log.(v_opt) atol=2e-3
        end
    end

    # @testset "Logistic regression via geometric programming" begin
    #     #=
    #     Section 7.1.1, Boyd & Vandenberghe, "Convex Optimization"
    #     cvxr.com/cvx/examples/
    #
    #     Solves the logistic regression problem re-formulated as a GP.
    #     The original log regression problem is:
    #
    #     minimize   sum_i(theta'*x_i) + sum_i( log(1 + exp(-theta'*x_i)) )
    #
    #     where x are explanatory variables and theta are model parameters.
    #     The equivalent GP is obtained by the following change of variables:
    #     z_i = exp(theta_i). The log regression problem is then a GP:
    #
    #     minimize   prod( prod(z_j^x_j) )*(prod( 1 + prod(z_j^(-x_j)) ))
    #
    #     with variables z and data x (explanatory variables).
    #     =#
    #
    #     a =  1
    #     b = -5
    #     M = 20
    #     u = 10*[0.7481934344558185,0.800187074451604,0.5526779186546091,0.7252363351360138,
    #         0.6610399174500099,0.41433761879760533,0.054666725368929026,0.205729269569759,
    #         0.930015247800057,0.5404894958921531,0.002089265140374774,0.11079906638513992,
    #         0.3011969126504366,0.3297993613893724,0.3208941434931467,0.42615481215185036,
    #         0.8872637536884693,0.953715521930721,0.8203465275422035,0.872935058440879]
    #     ybar = [0.8394192840775825,0.49753259430751684,0.687129915476693,0.6417040424559148,
    #         0.8792697149468107,0.5884540315572566,0.7235569366980186,0.008263405526384826,
    #         0.26010393902254947,0.23991239252986563,0.8999943204085568,0.010451240969007403,
    #         0.42624168523899475,0.9068096122478977,0.9545614303417385,0.8253383985140583,
    #         0.31854563622347576,0.24517545111977057,0.8227967366524322,0.4662118166849687]
    #     y = (ybar .< exp(a*u+b) ./ (1+exp(a*u+b)))
    #     ind_false = find(!y)
    #     ind_true  = find(y)
    #     # X is the sorted design matrix
    #     # first have true than false observations followed by the bias term
    #     X = vcat(u[ind_true], u[ind_false])
    #     X = hcat(X, ones(M))
    #     q = length(ind_true)
    #
    #     m = GPModel()
    #
    #     @variable(m, z[k=1:2])
    #     @variable(m, t[k=1:q])
    #     @variable(m, s[k=1:M])
    #
    #     @NLobjective(m, Min, prod{t[k], k=1:q}*prod{s[k], k=1:M})
    #
    #     @NLconstraint(m, lt[k=1:q], (z[1]^X[k,1])*(z[2]^X[k,2]) <= t[k])
    #     @NLconstraint(m, ls[k=1:M], 1 + (z[1]^(-X[k,1]))*(z[2]^(-X[k,2])) <= s[k])
    #
    #     # test continuous problem
    #     status = solve(m)
    #     @test status == :Optimal
    #     # comparing optimal continuous answers with cvx solutions
    #     @test getobjectivevalue(m) 2.150919077978715e4 1e-3
    #     @test log(getvalue(z[1])) 0.499270692481327 atol=1e-4
    #     @test log(getvalue(z[2])) -2.397533059296778 atol=1e-4
    # end

    @testset "Floor planning, $meth" for meth in methods
        #=
        Boyd, Kim, Vandenberghe, and Hassibi, "A Tutorial on Geometric Programming"

        Solves the problem of configuring and placing rectangles such
        that they do not overlap and that they minimize the area of the
        bounding box. This code solves the specific instances given
        in the GP tutorial. We have four rectangles with variable
        width w_i and height h_i. They need to satisfy area and aspect
        ration constraints. The GP is formulated as:

        minimize    max(wa+wb,wc+wd)*(max(ha,hb)+max(hc,hd))
            s.t.    wa*ha == area_a, wb*hb == area_b, ...
                    1/alpha_max <= ha/wa <= alpha_max, ...

        where variables are rectangle widths w's and heights h's.
        =#

        a = 0.2
        b = 0.5
        c = 1.5
        d = 0.5
        alpha = 2.11

        for solv in cont_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, wMax)
            @variable(m, habMax)
            @variable(m, hcdMax)
            @variable(m, wa)
            @variable(m, wb)
            @variable(m, wc)
            @variable(m, wd)
            @variable(m, ha)
            @variable(m, hb)
            @variable(m, hc)
            @variable(m, hd)

            @NLobjective(m, Min, wMax*(habMax + hcdMax))

            @NLconstraints m begin
                wa + wb <= wMax
                wc + wd <= wMax
                ha <= habMax
                hb <= habMax
                hc <= hcdMax
                hd <= hcdMax
                ha*wa == a
                hb*wb == b
                hc*wc == c
                hd*wd == d
                ha/wa <= alpha
                hb/wb <= alpha
                hc/wc <= alpha
                hd/wd <= alpha
                ha/wa >= 1/alpha
                hb/wb >= 1/alpha
                hc/wc >= 1/alpha
                hd/wd >= 1/alpha
            end

            # test continuous problem
            status = solve(m)

            @test status == :Optimal

            # comparing optimal answers with cvx solution
            @test getobjectivevalue(m) ≈ 2.929359634205730 atol=1e-4
            # @test getvalue(ha) 0.4067 atol=1e-4
            @test getvalue(hb) ≈ 0.4868 atol=1e-4
            @test getvalue(hc) ≈ 1.2247 atol=1e-4
            @test getvalue(hd) ≈ 1.0271 atol=1e-4
            # @test getvalue(wa) 0.4918 atol=1e-4
            @test getvalue(wb) ≈ 1.0271 atol=1e-4
            @test getvalue(wc) ≈ 1.2248 atol=1e-4
            @test getvalue(wd) ≈ 0.4868 atol=1e-4
        end
    end

    # @testset "Design of a cantilever beam" begin
    #     #=
    #     Boyd & Vandenberghe "Convex Optimization"
    #
    #     We have a segmented cantilever beam with N segments. Each segment
    #     has a unit length and variable width and height (rectangular profile).
    #     The goal is minimize the total volume of the beam, over all segment
    #     widths w_i and heights h_i, subject to constraints on aspect ratios,
    #     maximum allowable stress in the material, vertical deflection y, etc.
    #
    #         minimize    sum( w_i* h_i)
    #         s.t.        w_min <= w_i <= w_max,       for all i = 1,...,N
    #                     h_min <= h_i <= h_max
    #                     S_min <= h_i/w_i <= S_max
    #                     6*i*F/(w_i*h_i^2) <= sigma_max
    #                     6*F/(E*w_i*h_i^3) == d_i
    #                     (2*i - 1)*d_i + v_(i+1) <= v_i
    #                     (i - 1/3)*d_i + v_(i+1) + y_(i+1) <= y_i
    #                     y_1 <= y_max
    #
    #     with variables w_i, h_i, d_i, (i = 1,...,N) and v_i, y_i (i = 1,...,N+1).
    #     (Consult the book for other definitions and a recursive formulation of this problem.)
    #     =#
    #
    #     N = 8
    #     wmin = .1
    #     wmax = 100
    #     hmin = .1
    #     hmax = 6
    #     Smin = 1/5
    #     Smax = 5
    #     sigma_max = 1
    #     ymax = 10
    #     E = 1
    #     F = 1
    #
    #     m = GPModel()
    #
    #     @variable(m, wmin <= w[i=1:N] <= wmax)
    #     @variable(m, hmin <= h[i=1:N] <= hmax)
    #     @variable(m, v[i=1:N+1])
    #     @variable(m, y[i=1:N+1])
    #
    #     @NLobjective(m, Min, sum{w[i]*h[i], i=1:N})
    #
    #     @NLexpression(m, d[i=1:N], 6*F/(E*w[i]*h[i]^3))
    #     @NLconstraint(m, vle[i=1:N], (2*i-1)*d[i] + v[i+1] <= v[i])
    #     @NLconstraint(m, yle[i=1:N], (i-1/3)*d[i] + v[i+1] + y[i+1] <= y[i])
    #
    #     @NLconstraint(m, sMax[i=1:N], h[i]/w[i] <= Smax)
    #     @NLconstraint(m, sMin[i=1:N], h[i]/w[i] >= Smin)
    #     @NLconstraint(m, sig[i=1:N], 6*i*F/(w[i]*h[i]^2) <= sigma_max)
    #     @NLconstraint(m, yMax, y[1] <= ymax)
    #
    #     # test continuous problem
    #     status = solve(m)
    #     @test status == :Optimal
    #     # comparing optimal answers with cvx solution
    #     @test getobjectivevalue(m) 42.396550031186266 1e-5
    #     w_opt = [0.62144650818292,0.782973540379044,0.905967696541885,1.01240143291693,
    #         1.1003732685286,1.17622748621389,1.20000000091955,1.33333334198018]
    #     h_opt = [3.10723250879184,3.91486767920148,4.52983846322413,5.06200714644066,
    #         5.5018663252893,5.88113741411835,5.99999998429614,5.99999999183999]
    #     @test getvalue(w) w_opt 2e-5
    #     @test getvalue(h) h_opt 2e-5
    # end

    @testset "Frobenius norm diagonal scaling, $meth" for meth in methods
        #=
        Boyd & Vandenberghe "Convex Optimization"
        web.cvxr.com/cvx/examples/cvxbook/Ch04_cvx_opt_probs/html/frob_norm_diag_scaling.html

        Given a square matrix M, the goal is to find a vector (with dii > 0)
        such that ||DMD^{-1}||_F is minimized, where D = diag(d).
        The problem can be cast as an unconstrained geometric program:
            minimize sqrt( sum_{i,j=1}^{n} Mij^2*di^2/dj^2 )
        =#

        N = 4
        M = [0.6136879688613316 1.0335602629603533 -0.2950480421627936 0.5541491673222215
            -0.3978491302590482 1.178562648014158 -0.06526653255372962 0.21857690164142826
            -0.6997983901300008 -0.8699116880973898 1.1225766187758226 -2.654087224755067
            -0.40051651256252796 0.38234147083401476 -0.4742359774002496 -0.013651221368362231]

        for solv in cont_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, s)
            @variable(m, d[i=1:N])

            @NLobjective(m, Min, s^(1/2))

            @NLconstraint(m, s >= sum(M[i,j]^2*d[i]^2/d[j]^2 for i=1:N for j=1:N))

            # test continuous problem
            status = solve(m)

            @test status == :Optimal
            # comparing optimal answers with cvx solution
            @test getobjectivevalue(m) ≈ 2.746582079307096 atol=atol=1e-6
            # multiple opt solutions
            # d_opt = [0.861200327236396, 1.43343369595786, 0.476975798845936, 1.10877924241676]
            # @test getvalue(d) d_opt atol=1e-6
        end
    end

    @testset "Minimum spectral radius via Peron-Frobenius, $meth" for meth in methods
        #=
        Boyd & Vandenberghe "Convex Optimization"
        web.cvxr.com/cvx/examples/cvxbook/Ch04_cvx_opt_probs/html/min_spec_rad_ppl_dynamics.html

        The goal is to minimize the spectral radius of a square matrix A
        which is elementwise nonnegative, Aij >= 0 for all i,j. In this
        case A has a positive real eigenvalue lambda_pf (the Perron-Frobenius
        eigenvalue) which is equal to the spectral radius, and thus gives
        the fastest decay rate or slowest growth rate.

        We consider a specific example in which we want to find the fastest
        decay or slowest growth rate for the bacteria population governed
        by a simple dynamic model (see page 166). The problem is a GP:
            minimize   lambda
                s.t.   b1*v1 + b2*v2 + b3*v3 + b4*v4 <= lambda*v1
                      s1*v1 <= lambda*v2
                       s2*v2 <= lambda*v3
                       s3*v3 <= lambda*v4
                       1/2 <= ci <= 2
                       bi == bi^{nom}*(c1/c1^{nom})^alpha_i*(c2/c2^{nom})^beta_i
                       si == si^{nom}*(c1/c1^{nom})^gamma_i*(c2/c2^{nom})^delta_i

        with variables bi, si, ci, vi, lambda.
        =#

        c_nom = [1, 1]
        b_nom = [2, 3, 2, 1]
        alpha = [1, 1, 1, 1]
        beta  = [1, 1, 1, 1]
        s_nom = [1, 1, 3]
        gamma = [1, 1, 1]
        delta = [1, 1, 1]

        for solv in cont_solvers[meth]
            m = GPModel(method=meth, solver=solv)

            @variable(m, lambda)
            @variable(m, b[1:4])
            @variable(m, s[1:3])
            @variable(m, v[1:4])
            @variable(m, c[1:2])

            @NLobjective(m, Min, lambda)

            @NLconstraint(m, sum(b[i]*v[i] for i=1:4) <= lambda*v[1])
            @NLconstraint(m, s[1]*v[1] <= lambda*v[2])
            @NLconstraint(m, s[2]*v[2] <= lambda*v[3])
            @NLconstraint(m, s[3]*v[3] <= lambda*v[4])
            @NLconstraint(m, c[1] <= 2)
            @NLconstraint(m, c[2] <= 2)
            @NLconstraint(m, c[1] >= 0.5)
            @NLconstraint(m, c[2] >= 0.5)
            @NLconstraint(m, bn[i=1:4], b[i] == b_nom[i]*(c[1]/c_nom[1])^alpha[i]*(c[2]/c_nom[2])^beta[i])
            @NLconstraint(m, sn[i=1:3], s[i] == s_nom[i]*(c[1]/c_nom[1])^gamma[i]*(c[2]/c_nom[2])^delta[i])

            # test continuous problem
            status = solve(m)

            @test status == :Optimal
            # comparing optimal answers with cvx solution
            @test getobjectivevalue(m) ≈ 0.804067384844554 atol=1e-6
            b_opt = [0.500000000003113, 0.750000000004669, 0.500000000003113, 0.250000000001556]
            @test getvalue(b) ≈ b_opt atol=1e-6
            s_opt = [0.250000000001556, 0.250000000001556, 0.750000000004669]
            @test getvalue(s) ≈ s_opt atol=1e-6
            # multiple optimal v's
            # v_opt = [1.24367686858908, 0.386683035543248, 0.120227185822489, 0.112142826685870]
            # @test getvalue(v) v_opt atol=1e-6
            c_opt = [0.500000000001556, 0.500000000001556]
            @test getvalue(c) ≈ c_opt atol=1e-6
        end
    end
end
