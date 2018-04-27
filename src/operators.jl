#  Copyright 2016, Miles Lubin and contributors
#  Copyright 2018, Chris Coey and contributors
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
(+)(mon::Monomial) = 1*mon
(-)(mon::Monomial) = -1*mon
# --Number
(+)(mon::Monomial, x::Number) = mon + Monomial(x)
(-)(mon::Monomial, x::Number) = mon - Monomial(x)
(*)(mon::Monomial, x::Number) = Monomial(mon.c*x, mon.terms)
(/)(mon::Monomial, x::Number) = Monomial(mon.c/x, mon.terms)
(^)(m::Monomial, num::Integer) = Monomial(m.c^num, Dict{Int,Float64}([i => m.terms[i]*num for i in keys(m.terms)])) # for ambiguity
(^)(m::Monomial, num::Number) = Monomial(m.c^num, Dict{Int,Float64}([i => m.terms[i]*num for i in keys(m.terms)]))
# --Monomial
(+)(m_1::Monomial, m_2::Monomial) = Posynomial([m_1, m_2])
(-)(m_1::Monomial, m_2::Monomial) = Posynomial([m_1, -m_2])
function (*)(m_1::Monomial, m_2::Monomial)
    d = copy(m_1.terms)
    for (i,v) in m_2.terms
        d[i] = get(d, i, 0.0) + v
    end
    return Monomial(m_1.c*m_2.c, d)
end
function (/)(m_1::Monomial, m_2::Monomial)
    d = copy(m_1.terms)
    for (i,v) in m_2.terms
        d[i] = get(d, i, 0.0) - v
    end
    return Monomial(m_1.c/m_2.c, d)
end
# -- Posynomial
(+)(mon::Monomial,pos::Posynomial) = Posynomial(vcat(mon, pos.mons))
(-)(mon::Monomial,pos::Posynomial) = Posynomial(vcat(mon,map(-, pos.mons)))
(*)(mon::Monomial,pos::Posynomial) = Posynomial([mon*pm for pm in pos.mons])
(/)(mon::Monomial,pos::Posynomial) = error("Can't divide a Monomial by a Posynomial")

# Posynomial--
(+)(pos::Posynomial) = 1*pos
(-)(pos::Posynomial) = -1*pos
# --Number
(+)(pos::Posynomial, x::Number) = Posynomial(vcat(pos.mons, Monomial(x)))
(-)(pos::Posynomial, x::Number) = Posynomial(vcat(pos.mons, Monomial(-x)))
(*)(pos::Posynomial, x::Number) = Posynomial([pm*x for pm in pos.mons])
(/)(pos::Posynomial, x::Number) = Posynomial([pm/x for pm in pos.mons])
# --Monomial
(+)(pos::Posynomial, mon::Monomial) = Posynomial(vcat(pos.mons,mon))
(-)(pos::Posynomial, mon::Monomial) = Posynomial(vcat(pos.mons, -mon))
(*)(pos::Posynomial, mon::Monomial) = Posynomial([pm*mon for pm in pos.mons])
(/)(pos::Posynomial, mon::Monomial) = Posynomial([pm/mon for pm in pos.mons])
# --Posynomial
(+)(p_1::Posynomial, p_2::Posynomial) = Posynomial(vcat(p_1.mons, p_2.mons))
(-)(p_1::Posynomial, p_2::Posynomial) = Posynomial(vcat(p_1.mons, (-p_2).mons))
function (*)(p_1::Posynomial, p_2::Posynomial)
    mon_out = Monomial[]
    for m_1 in p_1.mons, m_2 in p_2.mons
        push!(mon_out, m_1*m_2)
    end
    return Posynomial(mon_out)
end
(/)(p_1::Posynomial, p_2::Posynomial) = error("Can't divide a Posynomial by a Posynomial")
