function combs(n,s)
    all_Jz=Vector{Any}()    
    function recurse(start,subset)
        if length(subset)==s
            Jz=subset;
            push!(all_Jz,Jz)
        end
        for j=start:n
            recurse(j+1,[subset...,j])
        end
    end
    recurse(1,[])
    return all_Jz
end
