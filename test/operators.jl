#  Copyright 2016, Miles Lubin and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using JuGP, JuMP
import JuGP: Monomial, Posynomial
if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

@testset "Operators" begin
@testset "Number--" begin
    @testset "--Monomial" begin
        @test 2 + Monomial(1) ==
            Posynomial([Monomial(2),Monomial(1)])
        @test 2 - Monomial(1) ==
            Posynomial([Monomial(2),Monomial(-1)])
        @test 2 * Monomial(3) == Monomial(6)
        @test 3 / Monomial(3,Dict(2=>2)) ==
            Monomial(1,Dict(2=>-2))
    end
    @testset "--Posynomial" begin
        @test 3 + (2 + Monomial(1)) ==
            Posynomial([Monomial(3),Monomial(2),Monomial(1)])
        @test 3 - (2 + Monomial(1)) ==
            Posynomial([Monomial(3),Monomial(-2),Monomial(-1)])
        @test 3 * (2 + Monomial(1)) ==
            Posynomial([Monomial(6),Monomial(3)])
        @test_throws ErrorException 3/(2 + Monomial(1))
    end
end
end
