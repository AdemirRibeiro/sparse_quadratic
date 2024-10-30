function gap(num_prob,tol,nfig)
    sv = true; szft=10;
    cor0 = RGB(.25,.25,.25); cor1 = RGB(.9,.1,.1); cor2 = RGB(.1,.1,.9);
    cor3 = RGB(.3,.9,.3); cor4 = RGB(.9,.1,.9); cor5 = RGB(.2,.5,.1);
    all_ρ = Vector{Any}();
    all_d_opt = Vector{Any}();
    all_p_opt = Vector{Any}();
    all_ρ = Vector{Any}();
    all_N = Vector{Any}();
    output = [];
    k = 0;
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
        N1 = factorial(n,n-s)/factorial(s);
        N = Int(N1);
        push!(all_N,N)
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
                    go = false
                    p_opt = pz;
                    Jz_opt = Jz;
                    ρ = executed_comb/N;
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
        (p_opt,indexmin) = findmin(all_pz);
        push!(all_p_opt,p_opt)
    end
    index_less = findall(all_ρ.<1);
    np = 1;
    aux = (plot([np;np],[all_d_opt[np];all_p_opt[np]],color=cor0,
    line=:dot,legend=false,xtickfont=font(szft),ytickfont=font(szft)))
    for np=2:num_prob
        aux = (plot!([np;np],[all_d_opt[np];all_p_opt[np]],color=cor0,
        line=:dot))
    end
    aux = (scatter!([1:num_prob],[all_d_opt],color=cor2,markershape=:star5);
           scatter!([1:num_prob],[all_p_opt],color=cor1,
                    xlims=(.5,num_prob+.5),legend=false);
           scatter!([index_less],[all_p_opt[index_less]],
           color=cor3);)
    display(aux)
    if sv
        nf = nfig; savefig("gap_values$nf.pdf");
    end
    aux = (scatter([1:num_prob],[all_ρ],color=cor0,legend=false,
                   xlims=(.5,num_prob+.5),ylims=(-.15,1.15),
                   xtickfont=font(szft),ytickfont=font(szft));
           scatter!([index_less],[all_ρ[index_less]],color=cor3))
    display(aux)
    if sv
        nf = nfig; savefig("gap_zero_proportion$nf.pdf");
    end
    ηs = round.(all_η, sigdigits=5);
    ps = round.(all_p_opt, sigdigits=5);
    ds = round.(all_d_opt, sigdigits=5);
    ρs = round.(all_ρ, sigdigits=5);
    output = [1:num_prob all_n all_s all_m ηs ps ds all_N ρs];
    return output
end
