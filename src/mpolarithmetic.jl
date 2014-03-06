function +(ps::MPol...)
    pu = PolUnion(ps...)
    r = zero(pu)
    for pindex = 1:length(ps)
        for (m, c) in ps[pindex]
            r[exps(pu, m, pindex)] += c
        end
    end
    r
end
function +{T<:Number,L<:Number}(p::MPol{T}, s::L)
    r = convert(MPol{promote_type(T, L)}, p)
    r[zeros(Int, nvars(p))] += s
    r
end
+{T<:Number}(s::T, p::MPol) = p + s


function -{T<:Number,L<:Number}(p::MPol{T}, q::MPol{L})
    pu = PolUnion(p, q)
    r = zero(pu)
    for (m, c) in p
        r[exps(pu, m, 1)] += c
    end
    for (m, c) in q
        r[exps(pu, m, 2)] -= c
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
    for (m, c) in p
        r[m] = -c
    end
    r
end
-{T<:Number}(s::T, p::MPol) = -p + s


function *{T<:Number,U<:Number}(p1::MPol{T}, p2::MPol{U})
    pu = PolUnion(p1, p2)
    r = zero(pu)
    for (m1, c1) in p1
        for (m2, c2) in p2
            r[exps(pu, m1, 1) + exps(pu, m2, 2)] += c1 * c2
        end
    end
    r
end
function *{T<:Number,U<:Number}(s::T, p::MPol{U})
    r = zero(MPol{promote_type(T,U)}, vars(p))
    for (m, c) in p
        r[m] = s * c
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
    for (m, c) in p
        r[m] = c/s
    end
    r
end


function conj(p::MPol)
    r = zero(p)
    for (m, c) in p
        r[m] = conj(c)
    end
    r
end