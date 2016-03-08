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
@testset "Monomial--" begin
    @testset "unary" begin
        @test +Monomial(1) == Monomial( 1)
        @test -Monomial(1) == Monomial(-1)
    end
    @testset "--Number" begin
        @test Monomial(1) + 2 ==
            Posynomial([Monomial(1),Monomial(2)])
        @test Monomial(1) - 1 ==
            Posynomial([Monomial(1),Monomial(-1)])
        @test Monomial(3) * 2 == Monomial(6)
        @test Monomial(3,Dict(2=>2)) / 3   == Monomial(1,Dict(2=>2))
        @test Monomial(3,Dict(2=>2)) ^ 2   == Monomial(9,Dict(2=>4))
        @test Monomial(3,Dict(2=>2)) ^ 2.0 == Monomial(9,Dict(2=>4))
    end
    @testset "--Monomial" begin
        @test Monomial(1) + Monomial(2) ==
            Posynomial([Monomial(1),Monomial(2)])
        @test Monomial(1) - Monomial(2) ==
            Posynomial([Monomial(1),Monomial(-2)])
        @test Monomial(2,Dict(1=>2,2=>3)) *
                Monomial(3,Dict(2=>4,3=>5)) ==
            Monomial(6,Dict(1=>2,2=>7,3=>5))
        @test Monomial(2,Dict(1=>2,2=>3)) /
                Monomial(3,Dict(2=>4,3=>5)) ==
            Monomial(2/3,Dict(1=>2,2=>-1,3=>-5))
    end
end
end
