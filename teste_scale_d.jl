function teste_scale_d(num_prob,nfig)
    sv = true;
    cor0 = RGB(.25,.25,.25); cor1 = RGB(.9,.1,.1); cor2 = RGB(.1,.1,.9);
    cor3 = RGB(.3,.9,.3); cor4 = RGB(.9,.1,.9); cor5 = RGB(.1,.5,.1);
    all_d_opt = Vector{Any}(); # using prog_maxmin
    all_d_opt2 = Vector{Any}(); # using prog_maxmin2
    all_d_opt3 = Vector{Any}(); # using prog_maxmin3
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
        d_opt = prog_maxmin(n,s,m,η,c,V,λ,A,b);
        d_opt2 = prog_maxmin2(n,s,m,η,c,V,λ,A,b);
        d_opt3 = prog_maxmin3(n,s,m,η,c,V,λ,A,b);
        push!(all_d_opt,d_opt)
        push!(all_d_opt2,d_opt2)
        push!(all_d_opt3,d_opt3)
    end
    indices = sortperm(all_d_opt);
    all_d_opt_sorted = all_d_opt[indices];
    all_d_opt2_sorted = all_d_opt2[indices];
    all_d_opt3_sorted = all_d_opt3[indices];
    aux = (scatter([1:num_prob],[all_d_opt_sorted],color=cor2,xlims=(.5,num_prob+.5),legend=false);
           scatter!([1:num_prob],[all_d_opt2_sorted],color=cor4);
           scatter!([1:num_prob],[all_d_opt3_sorted],color=cor5)
           )
    display(aux)
    if sv
        nf = nfig; savefig("/home/ademir/Ademir/artigos/Mael/teste_scale_d$nf.pdf");
    end
    return aux_seed'#, all_d_opt, all_d_opt2, all_d_opt3
end
