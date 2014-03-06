function printmonomial(io::IO, exps, vars)
    for i = 1:length(exps)        
        if exps[i] >= 1
            print(io, ' ')        
            print(io, vars[i])
            if exps[i] >= 2
                print(io, '^')
                print(io, exps[i])
            end
        end
    end
end

function show{T<:Number}(io::IO, p::MPol{T})
    first = true
    for (exps, c) in p
        if abs(c) > 0.00000001
        if first
            if typeof(c) <: Complex
                print(io, "($c)")
            elseif typeof(c) <: Rational
                print(io, "$(c)")
            else
                print(io, "$c")
            end
            printmonomial(io, exps, vars(p))
            
            first = false
        else
            if typeof(c) <: Complex
                print(io, " + ($c)")
            else
                if c >= zero(T)
                    print(io, " + $c")
                else
                    print(io, " - $(-c)")
                end
               # if typeof(c) <: Rational
               #     print(io, "*")
               # end
            end
            printmonomial(io, exps, vars(p))
        end
        end
    end
    if first
        print(io, zero(T))
    end
end