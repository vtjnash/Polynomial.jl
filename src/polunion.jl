immutable PolUnion
    newvars::Vector{Var}
    m::Array{Int, 2}
    eltype::Type
end

function PolUnion(ps::MPol...)
    newvars = isempty(ps) ? Var[] : union([vars(p) for p in ps]...)
    m = Array(Int, length(newvars), length(ps))
    for pindex = 1:length(ps)
        for i = 1:length(newvars)
            m[i, pindex] = findfirst(vars(ps[pindex]), newvars[i])
        end
    end
    PolUnion(newvars, m, promote_type([eltype(p) for p in ps]...))
end

eltype(pu::PolUnion) = pu.eltype

newvars(pu::PolUnion) = pu.newvars

nnewvars(pu::PolUnion) = length(newvars(pu))

newexps(pu::PolUnion, exps::Vector{Int}, polindex::Int) =
    [pu.m[i, polindex] == 0 ? 0 : exps[pu.m[i, polindex]] for i = 1:nnewvars(pu)]
    
newtype(pu::PolUnion) = MPol{eltype(pu)}

