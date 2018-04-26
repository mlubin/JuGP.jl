#  Copyright 2016, Miles Lubin and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

function check_expr_gp(ex::Expr)
    # println("Check: ", ex)

    # Process expression
    function descend(ex)
        if ex.head == :ref
            # A variable - simplest monomial
            mon = Monomial(1.0, Dict{Int,Float64}(ex.args[2] => 1.0))
            # println("Mon: ", mon)
            return mon
        else
            # First arg is operation type
            op = ex.args[1]
            # Collect any Xials (or numbers)
            xials = Union{Xial,Number}[]
            for a in ex.args
                if isa(a, Expr)
                    push!(xials, descend(a))
                elseif isa(a, Number)
                    push!(xials, a)
                end
            end
            # Merge them
            merged = reduce(eval(op), xials)
            # println("Merged: ", merged)
            return merged
        end
    end

    final = descend(ex)
    # println("Final: ", final)
    # println("")
    return final
end

function extract_constants(pos::Posynomial)
    keep_mons = Monomial[]
    constant = 0.0
    for mon in pos.mons
        if length(mon.terms) == 0
            constant += mon.c
        else
            push!(keep_mons, mon)
        end
    end
    return (Posynomial(keep_mons), constant)
end

extract_constants(mon::Monomial) = (mon, 0.0)
