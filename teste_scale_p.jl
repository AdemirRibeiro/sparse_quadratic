function teste_scale_p(num_prob,nfig)
    sv = true;
    cor0 = RGB(.25,.25,.25); cor1 = RGB(.9,.1,.1); cor2 = RGB(.1,.1,.9);
    cor3 = RGB(.3,.9,.3); cor4 = RGB(.9,.1,.9); cor5 = RGB(.1,.5,.1);
    all_pz = Vector{Any}();
    all_pz2 = Vector{Any}();
    all_pz3 = Vector{Any}();
    λ = 0;
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
        total_comb = factorial(n,n-s)/factorial(s);
        i = Int(rand(1:total_comb));
        Jz = combs(n,s)[i];
        pz = quadprog_max(n,s,m,η,c,V,λ,A,b,Jz)
        pz2 = quadprog_max2(n,s,m,η,c,V,λ,A,b,Jz)
        pz3 = quadprog_max3(n,s,m,η,c,V,λ,A,b,Jz)
        push!(all_pz,pz)
        push!(all_pz2,pz2)
        push!(all_pz3,pz3)
    end
    indices = sortperm(all_pz);
    all_pz_sorted = all_pz[indices];
    all_pz2_sorted = all_pz2[indices];
    all_pz3_sorted = all_pz3[indices];
    aux = (scatter([1:num_prob],[all_pz_sorted],color=cor2,xlims=(.5,num_prob+.5),legend=false);
           # scatter!([1:num_prob],[all_pz2_sorted],color=cor4);
           scatter!([1:num_prob],[all_pz3_sorted],color=cor5)
           )
    display(aux)
    if sv
        nf = nfig; savefig("/home/ademir/Ademir/artigos/Mael/teste_scale_p$nf.pdf");
    end
    return all_pz, all_pz2, all_pz3
end
