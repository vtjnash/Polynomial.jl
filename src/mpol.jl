immutable MPol{T<:Number}
    coeffs::Dict{Vector{Int},T}
    vars::Vector{Var}
end

MPol{T<:Number}(coeffs::Dict{Vector{Int},T}, vars::Vector{Var}) = MPol{T}(coefss, vars)
MPol{T<:Number}(::Type{T}, vars::Vector{Var}) = MPol{T}(Dict{Vector{Int},T}(), vars)

eltype{T<:Number}(::MPol{T}) = T

zero{T<:Number}(::Type{MPol{T}}) = MPol(T, Var[])
zero{T<:Number}(::Type{MPol{T}}, vars::Vector{Var}) = MPol(T, vars)
zero{T<:Number}(p::MPol{T}) = zero(MPol{T}, vars(p))

one{T<:Number}(::Type{MPol{T}}) = MPol{T}([Int[] => one(T)], Var[])
one{T<:Number}(::Type{MPol{T}}, vars::Vector{Var}) = MPol{T}([zeros(Int, length(vars)) => one(T)], vars)
one{T<:Number}(p::MPol{T}) = one(MPol{T}, vars(p))

monomial{T<:Number}(::Type{MPol{T}}, var::Var) = MPol{T}([[1] => one(T)], [var])

vars(p::MPol) = p.vars
vars(p::MPol...) = newvars(PolUnion(p...))

nvars(p::MPol) = length(vars(p))

coeffs(p::MPol) = p.coeffs

promote_rule{T<:Number,L<:Number}(::Type{MPol{T}}, ::Type{MPol{L}}) = MPol{promote_type(T, L)}
promote_rule{T<:Number,L<:Number}(::Type{MPol{T}}, ::Type{L}) = MPol{promote_type(T, L)}

convert{T<:Number,L<:Number}(::Type{MPol{T}}, c::L) = MPol{T}([Int[] => convert(T, c)], []) 
function convert{T<:Number,L<:Number}(::Type{MPol{T}}, p::MPol{L})
    r = MPol(T, vars(p))
    for (m, c) in p
        r[m] = convert(T, c)
    end
    r
end

monomials(p::MPol) = keys(coeffs(p))

getindex{T<:Number}(p::MPol{T}, exps::Vector{Int}) = 
    get(coeffs(p), exps, zero(T))

getindex{T<:Number}(p::MPol{T}, exps::Int...) = p[[exps...]]

function setindex!{T<:Number}(p::MPol{T}, v::T, m::Int...)
    if isapprox(v, zero(T))
        delete!(coeffs(p), [m...])
    else
        coeffs(p)[[m...]] = v
    end
end

setindex!{T<:Number}(p::MPol{T}, v::T, m::Vector{Int}) =
    p[m...] = v
    
start(p::MPol) = start(coeffs(p))
next(p::MPol, state) = next(coeffs(p), state)
done(p::MPol, state) = done(coeffs(p), state)

function copy{T<:Number}(p::MPol{T})
    r = MPol(T, vars(p))
    for (m, c) in p
        r[m] = c
    end
    r
end

function deg(p::MPol)
    d = 0
    for exps in monomials(p)
        s = sum(exps)
        if s > d
            d = s
        end
    end
    d
end

function ==(p1::MPol, p2::MPol)
    for (m, c) in p1
        if p2[m] != c
            return false
        end
    end
    for (m, c) in p2
        if p1[m] != c
            return false
        end
    end
    true
end