#  Copyright 2016, Miles Lubin and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Number--
# --Monomial
(+)(x::Number, mon::Monomial) = Monomial(x) + mon
(-)(x::Number, mon::Monomial) = Monomial(x) - mon
(*)(x::Number, mon::Monomial) = Monomial(x*mon.c, mon.terms)
(/)(x::Number, mon::Monomial) = x*mon^-1
# --Posynomial
(+)(x::Number, pos::Posynomial) = Monomial(x) + pos
(-)(x::Number, pos::Posynomial) = Monomial(x) - pos
(*)(x::Number, pos::Posynomial) = Posynomial(map(mon->(x*mon), pos.mons))
(/)(x::Number, pos::Posynomial) = error("Can't divide a scalar by a Posynomial")


# Mon-number
(*)(mon::Monomial,num::Number) = num*mon
function (/)(num::Number, m::Monomial)
    return Monomial(num/m.c, Dict{Int,Float64}(map(x->(x[1],-x[2]),m.terms)))
end
(-)(m::Monomial, num::Number) = m - Monomial(num)
(+)(m::Monomial, num::Number) = m + Monomial(num)
# for ambiguity
(^)(m::Monomial, num::Integer) = Monomial(m.c^num,
                                    Dict{Int,Float64}([i => m.terms[i]*num for i in keys(m.terms)]))
(^)(m::Monomial, num::Number) = Monomial(m.c^num,
                                    Dict{Int,Float64}([i => m.terms[i]*num for i in keys(m.terms)]))

# Mon-Mon
(+)(m::Monomial, n::Monomial) = Posynomial([m,n])
(-)(m::Monomial, n::Monomial) = Posynomial([m,-1*n])
function (*)(m::Monomial, n::Monomial)
    d = copy(m.terms)
    for (i,v) in n.terms
        d[i] = get(d,i,0.0) + v
    end
    return Monomial(m.c*n.c, d)
end
function (/)(m::Monomial, n::Monomial)
    d = copy(m.terms)
    for (i,v) in n.terms
        d[i] = get(d,i,0.0) - v
    end
    return Monomial(m.c/n.c, d)
end


(-)(pos::Posynomial, num::Number) = pos - Monomial(num)
(-)(pos::Posynomial, mon::Monomial) = pos + (-1)*mon
(-)(mon::Monomial, pos::Posynomial) = mon + (-1)*pos
function (*)(pos::Posynomial, m::Monomial)
    return Posynomial([mon*m for mon in pos.mons])
end
(*)(mon::Monomial,pos::Posynomial) = pos*mon

(+)(p1::Posynomial, m::Monomial) = Posynomial(vcat(p1.mons,m))
(+)(mon::Monomial,pos::Posynomial) = Posynomial(vcat(mon,pos.mons))

(+)(p1::Posynomial, p2::Posynomial) = Posynomial([p1.mons;p2.mons])
