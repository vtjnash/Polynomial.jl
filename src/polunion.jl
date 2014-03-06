immutable PolUnion
    vars::Vector{Var}
    mat::Array{Int, 2}
    eltype::Type
end

function PolUnion(ps::MPol...)
    newvars = isempty(ps) ? Var[] : union([vars(p) for p in ps]...)
    mat = Array(Int, length(newvars), length(ps))
    for pindex = 1:length(ps)
        for i = 1:length(newvars)
            mat[i, pindex] = findfirst(vars(ps[pindex]), newvars[i])
        end
    end
    PolUnion(newvars, mat, promote_type([eltype(p) for p in ps]...))
end

eltype(pu::PolUnion) = pu.eltype

vars(pu::PolUnion) = pu.vars

nvars(pu::PolUnion) = length(vars(pu))

exps(pu::PolUnion, exps::Vector{Int}, polindex::Int) =
    [pu.mat[i, polindex] == 0 ? 0 : exps[pu.mat[i, polindex]] for i = 1:nvars(pu)]
    
zero(pu::PolUnion) = zero(MPol{eltype(pu)}, vars(pu))

one(pu::PolUnion) = one(MPol{eltype(pu)}, vars(pu))