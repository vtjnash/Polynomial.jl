module Polynomial

export 
    coeffs,
    UPol, polyval, polyint, polyder, poly, roots, degree, var, monomial,
    MPol, vars

import Base: zero, one
import Base: show, print, length, endof, getindex, setindex!, copy, promote_rule, convert, start, next, done, eltype
import Base: *, /, //, -, +, ==, divrem, conj, rem

typealias Var Symbol

eps{T}(::Type{T}) = convert(T,0)
eps{F<:FloatingPoint}(x::Type{F}) = Base.eps(F)
eps{T}(x::Type{Complex{T}}) = eps(T)

include("upol.jl")

include("mpol.jl")
include("mpolarithmetic.jl")
include("mpolprinting.jl")
include("polunion.jl")

end # module