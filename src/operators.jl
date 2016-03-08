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

# Monomial--
(+)(mon::Monomial) =  1*mon
(-)(mon::Monomial) = -1*mon
# --Number
(+)(mon::Monomial, x::Number) = mon + Monomial(x)
(-)(mon::Monomial, x::Number) = mon - Monomial(x)
(*)(mon::Monomial, x::Number) = Monomial(mon.c*x, mon.terms)
(/)(mon::Monomial, x::Number) = Monomial(mon.c/x, mon.terms)
(^)(m::Monomial, num::Integer) = Monomial(m.c^num,  # for ambiguity
    Dict{Int,Float64}([i => m.terms[i]*num for i in keys(m.terms)]))
(^)(m::Monomial, num::Number) = Monomial(m.c^num,
    Dict{Int,Float64}([i => m.terms[i]*num for i in keys(m.terms)]))
# --Monomial
(+)(m_1::Monomial, m_2::Monomial) = Posynomial([m_1, m_2])
(-)(m_1::Monomial, m_2::Monomial) = Posynomial([m_1,-m_2])
function (*)(m_1::Monomial, m_2::Monomial)
    d = copy(m_1.terms)
    for (i,v) in m_2.terms
        d[i] = get(d,i,0.0) + v
    end
    return Monomial(m_1.c*m_2.c, d)
end
function (/)(m_1::Monomial, m_2::Monomial)
    d = copy(m_1.terms)
    for (i,v) in m_2.terms
        d[i] = get(d,i,0.0) - v
    end
    return Monomial(m_1.c/m_2.c, d)
end
# -- Posynomial
(+)(mon::Monomial,pos::Posynomial) = Posynomial(vcat(mon,pos.mons))
(-)(mon::Monomial,pos::Posynomial) = Posynomial(vcat(mon,map(-,pos.mons)))
(*)(mon::Monomial,pos::Posynomial) = Posynomial([mon*pm for pm in pos.mons])
(/)(mon::Monomial,pos::Posynomial) = error("Can't divide a Monomial by a Posynomial")


(-)(pos::Posynomial, num::Number) = pos - Monomial(num)
(-)(pos::Posynomial, mon::Monomial) = pos + (-1)*mon
(*)(pos::Posynomial, m::Monomial) =m * pos
(+)(p1::Posynomial, m::Monomial) = Posynomial(vcat(p1.mons,m))
(+)(p1::Posynomial, p2::Posynomial) = Posynomial([p1.mons;p2.mons])
