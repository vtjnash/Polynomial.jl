immutable UPol{T<:Number}
    a::Vector{T}
    nzfirst::Int #for effiencicy, track the first non-zero index
    var::Symbol
    function UPol(a::Vector{T}, var::Var=:x)
        nzfirst = 0 #find and chop leading zeros
        for i = 1:length(a)
            if abs(a[i]) > 2*eps(T)
                break
            end
            nzfirst = i
        end
        new(a, nzfirst, var)
    end
end

UPol{T}(a::Vector{T}, var::Var=:x) = UPol{T}(a, var)

convert{T}(::Type{UPol{T}}, p::UPol) = UPol(convert(Vector{T}, p.a), var(p))
promote_rule{T, S}(::Type{UPol{T}}, ::Type{UPol{S}}) = UPol{promote_type(T, S)}
eltype{T}(::UPol{T}) = T

length(p::UPol) = length(p.a)-p.nzfirst
endof(p::UPol) = length(p)
deg(p::UPol) = length(p) - 1

getindex(p::UPol, i) = p.a[i+p.nzfirst]
setindex!(p::UPol, v, i) = (p.a[i+p.nzfirst] = v)

copy(p::UPol) = UPol(copy(p.a[1+p.nzfirst:end]), var(p))

var(p::UPol) = p.var
vars(p::UPol) = [var(p)]

zero{T}(p::UPol{T}) = UPol([zero(T)], var(p))
zero{T}(::Type{UPol{T}}) = UPol([zero(T)])
zero{T}(::Type{UPol{T}}, var::Var) = UPol([zero(T)], var)

one{T}(p::UPol{T}) = UPol([one(T)], var(p))
one{T}(::Type{UPol{T}}) = UPol([one(T)])
one{T}(::Type{UPol{T}}, var::Var) = UPol([one(T)], var)

monomial{T}(p::UPol{T}) = UPol([one(T), zero(T)], var(p))
monomial{T}(::Type{UPol{T}}, var::Var=:x) = UPol([one(T), zero(T)], var)

function show(io::IO, p::UPol)
    print(io,"UPol(")
    print(io,p)
    print(io,")")
end

function print{T}(io::IO, p::UPol{T})
    n = length(p)
    if n <= 0
        print(io,"0")
    else
        for j = 1:n
            pj = p[j]
            magpj = abs(pj)
            if magpj > 2*eps(T)
                if j == 1 
                    pj < 0 && print(io, "-")    #Prepend - if first and negative
                else
                    pj < 0 ? print(io," - ") : print(" + ")
                end
                #Print pj if pj is the last coefficient, or pj is not identically 1
                if j == n || abs(magpj - 1) > 2*eps(T)
                    print(io, magpj)
                end
                exp = n-j
                if exp > 0
                    print(io, var(p))
                    if exp > 1
                        print(io, '^', exp)
                    end
                end
            end
        end
    end
end

function print{T<:Complex}(io::IO, p::UPol{T})
    n = length(p)
    if n <= 0
        print(io,"0")
    else
        for j = 1:n
            pj = p[j]
            abs_repj = abs(real(pj))
            abs_impj = abs(imag(pj))
            if abs(pj) > 2*eps(T)
                if !(abs_impj > 2*eps(T))
                    if j > 1 
                        real(pj) < 0 ? print(io," - ") : print(" + ")
                    else
                        real(pj) < 0 && print(io, "-")    #Prepend - if first and negative
                    end
                else
                    j > 1 && print(io, " + ")
                end
                if abs_repj > 2*eps(T)    #Real part is not 0
                    if abs_impj > 2*eps(T)    #Imag part is not 0
                        print(io,'(',pj,')')
                    elseif abs(abs_repj - 1) > 2*eps(T) || j == n 
                        print(io,abs_repj)
                    end
                else
                    if abs_impj > 2*eps(T)
                        print(io,'(', imag(pj),"im)")
                    end
                end
                exp = n-j
                if exp > 0
                    print(io, var(p))
                    if exp > 1
                        print(io, '^', exp)
                    end
                end
            end
        end
    end
end

*(c::Number, p::UPol) = UPol(c * p.a[1+p.nzfirst:end], var(p))
*(p::UPol, c::Number) = UPol(c * p.a[1+p.nzfirst:end], var(p))
/(p::UPol, c::Number) = UPol(p.a[1+p.nzfirst:end] / c, var(p))
-(p::UPol) = UPol(-p.a[1+p.nzfirst:end], var(p))

-(p::UPol, c::Number) = +(p, -c)
+(c::Number, p::UPol) = +(p, c)
function +(p::UPol, c::Number)
    if length(p) < 1
        return UPol([c,], var(p))
    else
        p2 = copy(p)
        p2.a[end] += c
        return p2
    end
end
function -(c::Number, p::UPol)
    if length(p) < 1
        return UPol([c,], var(p))
    else
        p2 = -p
        p2.a[end] += c
        return p2
    end
end

function +{T,S}(p1::UPol{T}, p2::UPol{S})
    if p1.var != p2.var
        error("Polynomials must have same variable")
    end
    R = promote_type(T,S)
    n = length(p1)
    m = length(p2)
    if n > m
        a = Array(R, n)
        for i = 1:m
            a[n-m+i] = p1[n-m+i] + p2[i]
        end
        for i = 1:n-m
            a[i] = p1[i]
        end
    else
        a = Array(R, m)
        for i = 1:n
            a[m-n+i] = p1[i] + p2[m-n+i]
        end
        for i = 1:m-n
            a[i] = p2[i]
        end
    end
    UPol(a, p1.var)
end

function -{T,S}(p1::UPol{T}, p2::UPol{S})
    if p1.var != p2.var
        error("Polynomials must have same variable")
    end
    R = promote_type(T,S)
    n = length(p1)
    m = length(p2)
    if n > m
        a = Array(R, n)
        for i = 1:m
            a[n-m+i] = p1[n-m+i] - p2[i]
        end
        for i = 1:n-m
            a[i] = p1[i]
        end
    else
        a = Array(R, m)
        for i = 1:n
            a[m-n+i] = p1[i] - p2[m-n+i]
        end
        for i = 1:m-n
            a[i] = -p2[i]
        end
    end
    UPol(a, p1.var)
end

function *{T,S}(p1::UPol{T}, p2::UPol{S})
    if p1.var != p2.var
        error("Polynomials must have same variable")
    end
    R = promote_type(T,S)
    n = length(p1)
    m = length(p2)
    if n == 0 || m == 0
        return UPol(R[], p1.var)
    end
    a = zeros(R, n+m-1)
    for i = 1:length(p1)
        for j = 1:length(p2)
            a[i+j-1] += p1[i] * p2[j]
        end
    end
    UPol(a, p1.var)
end

function divrem{T, S}(num::UPol{T}, den::UPol{S})
    if num.var != den.var
        error("UPolnomials must have same variable")
    end
    m = length(den)
    if m == 0
        throw(DivideError())
    end
    R = typeof(one(T)/one(S))
    n = length(num)
    deg = n-m+1
    if deg <= 0
        return zero(UPol{R}), convert(UPol{R}, num)
    end
    d = zeros(R, n)
    q = zeros(R, deg)
    r = zeros(R, n)
    r[:] = num.a[1+num.nzfirst:end]
    for i = 1:deg
        quot = r[i] / den[1]
        q[i] = quot
        if i > 1
            d[i-1] = 0
            r[i-1] = 0
        end
        for j = 1:m
            k = i+j-1
            elem = den[j]*quot
            d[k] = elem
            r[k] -= elem
        end
    end
    return UPol(q, num.var), UPol(r, num.var)
end
/(num::UPol, den::UPol) = divrem(num, den)[1]
rem(num::UPol, den::UPol) = divrem(num, den)[2]

function ==(p1::UPol, p2::UPol)
    if length(p1) != length(p2)
        return false
    elseif p1.var != p2.var
        return false
    else
        return p1.a[1+p1.nzfirst:end] == p2.a[1+p2.nzfirst:end]
    end
end

function polyval{T}(p::UPol{T}, x::Number)
    R = promote_type(T, typeof(x))
    lenp = length(p)
    if lenp == 0
        return zero(R)
    else
        y = convert(R, p[1])
        for i = 2:lenp
            y = p[i] + x.*y
        end
        return y
    end
end

polyval(p::UPol, v::AbstractVector) = map(x->polyval(p, x), v)

function polyint{T}(p::UPol{T}, k::Number=0)
    n = length(p)
    R = typeof(one(T)/1)
    a2 = Array(R, n+1)
    for i = 1:n
        a2[i] = p[i] / (n-i+1)
    end
    a2[end] = k
    UPol(a2, var(p))
end

function polyder{T}(p::UPol{T})
    n = length(p)
    if n > 0
        a2 = Array(T, n-1)
        for i = 1:n-1
            a2[i] = p[i] * (n-i)
        end
    else
        a2 = zeros(T, 0)
    end
    UPol(a2, var(p))
end

# create a UPol object from its roots
function poly{T}(r::AbstractVector{T}, var::Var=:x)
    n = length(r)
    c = zeros(T, n+1)
    c[1] = 1
    for j = 1:n
        c[2:j+1] = c[2:j+1]-r[j]*c[1:j]
    end
    return UPol(c, var)
end
poly(A::Matrix, var::Var=:x) = poly(eig(A)[1], var)

roots{T}(p::UPol{Rational{T}}) = roots(convert(UPol{promote_type(T, Float64)}, p))

# compute the roots of a polynomial
function roots{T}(p::UPol{T})
    R = promote_type(T, Float64)
    num_zeros = 0
    if length(p) == 0
        return zeros(R, 0)
    end
    while abs(p[end-num_zeros]) <= 2*eps(T)
        if num_zeros == length(p)-1
            return zeros(R, 0)
        end
        num_zeros += 1
    end
    n = length(p)-num_zeros-1
    if n < 1
        return zeros(R, length(p)-1)
    end
    companion = zeros(R, n, n)
    a0 = p[end-num_zeros]
    for i = 1:n-1
        companion[1,i] = -p[end-num_zeros-i] / a0
        companion[i+1,i] = 1;
    end
    companion[1,end] = -p[1] / a0
    D,V = eig(companion)
    r = zeros(eltype(D),length(p)-1)
    r[1:n] = 1./D
    return r
end