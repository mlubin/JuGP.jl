#  Copyright 2016, Miles Lubin and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# inspiration from Iain's initial work at https://github.com/IainNZ/GPTest

import Base: (*), (+), (-), (/), (^), (==)

abstract type Xial end

mutable struct Monomial <: Xial
    c::Float64
    terms::Dict{Int,Float64} # variable index to power (coefficient)
end

Monomial(c::Number) = Monomial(c, Dict{Int,Float64}())

Base.convert(::Type{Monomial}, c::Number) = Monomial(c)

(==)(m1::Monomial, m2::Monomial) = (m1.c == m2.c) && (m1.terms == m2.terms)

function Base.print(io::IO, mon::Monomial)
    print(io, mon.c)
    for (i,v) in mon.terms
        print(io, "*x_{$i}^{$v}")
    end
end

mutable struct Posynomial <: Xial
    mons::Vector{Monomial}
end

Posynomial(c::Number) = Posynomial(Monomial(c))
Posynomial(m::Monomial) = Posynomial(Monomial[m])
Posynomial(args...) = Posynomial(Monomial[args...])

function (==)(p1::Posynomial, p2::Posynomial)
    # Doesn't handle sorting differences
    if length(p1.mons) != length(p2.mons)
        return false
    end
    all(pair -> (pair[1] == pair[2]), zip(p1.mons, p2.mons))
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

mutable struct GPData
    discretevalues::Dict{Int,Vector{Float64}}
end

GPData() = GPData(Dict{Int,Vector{Float64}}())
