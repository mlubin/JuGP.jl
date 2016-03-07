#  Copyright 2016, Miles Lubin and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# inspiration from Iain's initial work at https://github.com/IainNZ/GPTest

import Base: (*), (+), (-), (/), (^)

abstract Xial

type Monomial <: Xial
    c::Float64
    terms::Dict{Int,Float64} # variable index to power (coefficient)
end
Monomial() = Monomial(1.0,Dict{Int,Float64}())
Monomial(c::Number) = Monomial(c, Dict{Int,Float64}())

function Base.print(io::IO, mon::Monomial)
    print(io, mon.c)
    for (i,v) in mon.terms
        print(io, "*x_{$i}^{$v}")
    end
end

type Posynomial <: Xial
    mons::Vector{Monomial}
end
function Base.print(io::IO, pos::Posynomial)    
    if length(pos.mons) == 0
        print(io, "0")
    elseif length(pos.mons) == 1
        print(io, pos.mons[1])
    else
        print(io, "[")
        for mon in pos.mons[1:end-1]
            print(io, mon, " + ")
        end
        print(io, pos.mons[end], "]")
    end
end

# Mon-number
function (*)(num::Number, m::Monomial)
    return Monomial(num*m.c, m.terms)
end
function (/)(num::Number, m::Monomial)
    return Monomial(num/m.c, Dict{Int,Float64}(map(x->(x[1],-x[2]),m.terms)))
end
(+)(m::Number, num::Monomial) = Monomial(num) + m
(-)(m::Monomial, num::Number) = m - Monomial(num)
(+)(m::Monomial, num::Number) = m + Monomial(num)
(-)(num::Number, m::Monomial) = Monomial(num) - m
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

(*)(num::Number, pos::Posynomial) = Posynomial(map(m->(num*m), pos.mons))
(-)(pos::Posynomial, num::Number) = pos - Monomial(num)
(-)(pos::Posynomial, m::Monomial) = pos + (-1)*m
function (*)(pos::Posynomial, m::Monomial)
    return Posynomial([mon*m for mon in pos.mons])
end

function (+)(p1::Posynomial, m::Monomial)
    return Posynomial([p1.mons;m])
end

function (+)(p1::Posynomial, p2::Posynomial)
    return Posynomial([p1.mons;p2.mons])
end
