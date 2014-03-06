function +(ps::MPol...)
    pu = PolUnion(ps...)
    r = zero(newtype(pu), newvars(pu))
    for pindex = 1:length(ps)
        for (exps, c) in ps[pindex]
            r[newexps(pu, exps, pindex)] += c
        end
    end
    r
end
function +{T<:Number,L<:Number}(p::MPol{T}, s::L)
    r = convert(MPol{promote_type(T, L)}, p)
    r[zeros(Int, nvars(p))] = r[zeros(Int, nvars(p))] + s
    r
end
+{T<:Number}(s::T, p::MPol) = p + s


function -{T<:Number,L<:Number}(p::MPol{T}, q::MPol{L})
    pu = PolUnion(p, q)
    r = zero(newtype(pu), newvars(pu))
    for (exps, c) in p
        r[newexps(pu, exps, 1)] += c
    end
    for (exps, c) in q
        r[newexps(pu, exps, 2)] -= c
    end
    r
end
function -{T<:Number,U<:Number}(p::MPol{T}, s::U)
    r = convert(MPol{promote_type(T, U)}, p)
    r[zeros(Int, nvars(p))] -= s
    r
end
function -{T<:Number}(p::MPol{T})
    r = MPol(T, vars(p))
    for (exps, c) in p
        r[exps] = -c
    end
    r
end
-{T<:Number}(s::T, p::MPol) = -p + s


function *{T<:Number,U<:Number}(p1::MPol{T}, p2::MPol{U})
    pu = PolUnion(p1, p2)
    r = zero(newtype(pu), newvars(pu))
    for (exps1, c1) in p1
        for (exps2, c2) in p2
            r[newexps(pu, exps1, 1) + newexps(pu, exps2, 2)] += c1 * c2
        end
    end
    r
end
function *{T<:Number,U<:Number}(s::T, p::MPol{U})
    r = zero(MPol{promote_type(T,U)}, vars(p))
    for (exps, c) in p
        r[exps] = s * c
    end
    r
end
*{T<:Number}(p::MPol, s::T) = s * p


function ^{T<:Number}(p::MPol{T}, power::Integer)
    @assert power >= 0
    if power == 0
        return one(p)
    elseif power == 1
        return p
    else
        f, r = divrem(power, 2)
        return p^(f+r) * p^f
    end
end


function /{T<:Number,U<:Number}(p::MPol{T}, s::U)
    r = zero(MPol{promote_type(T,U)}, vars(p))
    for (exps, c) in p
        r[exps] = c/s
    end
    r
end


function conj(p::MPol)
    r = zero(p)
    for (exps, c) in p
        r[exps] = conj(c)
    end
    r
end