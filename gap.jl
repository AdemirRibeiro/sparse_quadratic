function gap(num_prob,tol,nfig)
    sv = true; szft=10;
    cor0 = RGB(.25,.25,.25); cor1 = RGB(.9,.1,.1); cor2 = RGB(.1,.1,.9);
    cor3 = RGB(.3,.9,.3); cor4 = RGB(.9,.1,.9); cor5 = RGB(.2,.5,.1);
    all_ρ = Vector{Any}();
    all_d_opt = Vector{Any}(); # using prog_maxmin
    all_p_opt = Vector{Any}(); # using quadprog_max
    k = 0;
    all_ξ = Vector{Any}();
    d_opt = 0;
    Jz_opt = [];
    aux_seed, all_n, all_s, all_m, all_η, all_c, all_V, all_λ, all_A, all_b = data(num_prob);
    for np = 1:num_prob
        n = all_n[np];
        s = all_s[np];
        m = all_m[np];
        η = all_η[np];
        c = all_c[np];
        V = all_V[np];
        λ = all_λ[np];
        A = all_A[np];
        b = all_b[np];
        sqrtλ = sqrt.(λ);
        all_Jz=Vector{Any}();
        all_pz = Vector{Any}();
        total_comb=factorial(n,n-s)/factorial(s);
        executed_comb=0;
        d_opt = prog_maxmin(n,s,m,η,c,V,λ,A,b);
        push!(all_d_opt,d_opt)
        go = true
        function combination(start,subset)
            if length(subset) == s
                executed_comb = executed_comb + 1;
                Jz = subset;
                pz = quadprog_max(n,s,m,η,c,V,λ,A,b,Jz);
                push!(all_Jz,Jz)
                push!(all_pz,pz)
                if abs(pz-d_opt)/(abs(d_opt)) < tol
                    go=false
                    p_opt = pz;
                    Jz_opt = Jz;
                    ρ = executed_comb/total_comb;
                    push!(all_ρ,ρ)
                end
            end
            for j = start:n
                if go; combination(j+1,[subset...,j]); end
            end
        end
        combination(1,[])
        if go
            push!(all_ρ,1)
        end
        (p_opt,indexmin)=findmin(all_pz);
        push!(all_p_opt,p_opt)
    end
    indices = sortperm(all_p_opt);
    all_p_opt_sorted = all_p_opt[indices]; # The same (primal) optimal values, but sorted
    all_d_opt_sorted = all_d_opt[indices]; # The dual optimal values (not necessarily sorted)
    all_ρ_sorted = all_ρ[indices];
    index_less_sorted = findall(all_ρ_sorted.<1);
    np=1;
    aux = (plot([np;np],[all_d_opt_sorted[np];all_p_opt_sorted[np]],color=cor0,
    line=:dot,legend=false,xtickfont=font(szft),ytickfont=font(szft)))
    for np=2:num_prob
        aux = (plot!([np;np],[all_d_opt_sorted[np];all_p_opt_sorted[np]],color=cor0,
        line=:dot))
    end
    aux = (scatter!([1:num_prob],[all_d_opt_sorted],color=cor2,markershape=:star5);
           scatter!([1:num_prob],[all_p_opt_sorted],color=cor1,
                    xlims=(.5,num_prob+.5),legend=false);
           scatter!([index_less_sorted],[all_p_opt_sorted[index_less_sorted]],
           color=cor3);)
    display(aux)
    if sv
        nf = nfig; savefig("/home/ademir/Ademir/artigos/Mael/fig/gap_values_sorted$nf.pdf");
    end
    aux = (scatter([1:num_prob],[all_ρ_sorted],color=cor0,legend=false,
                   xlims=(.5,num_prob+.5),ylims=(-.15,1.15),
                   xtickfont=font(szft),ytickfont=font(szft));
           scatter!([index_less_sorted],[all_ρ_sorted[index_less_sorted]],color=cor3))
    display(aux)
    if sv
        nf = nfig; savefig("/home/ademir/Ademir/artigos/Mael/fig/gap_zero_proportion_sorted$nf.pdf");
    end
    return aux_seed, index_less_sorted
end
