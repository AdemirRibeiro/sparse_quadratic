function teste_quadprog(nfig)
    sv = true;
    cor0 = RGB(.25,.25,.25); cor1 = RGB(.9,.1,.1); cor2 = RGB(.1,.1,.9);
    cor3 = RGB(.3,.9,.3); cor4 = RGB(.9,.1,.9); cor5 = RGB(.1,.5,.1);
    all_pz = Vector{Any}();
    λ = 0;
    aux_seed, all_n, all_s, all_m, all_η, all_c, all_V, all_λ, all_A, all_b = data(1);
    np = 1;
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
    tc = Int(total_comb);
    all_Jz = combs(n,s);
    for i = 1:tc
        Jz = all_Jz[i];
        pz = quadprog_max(n,s,m,η,c,V,λ,A,b,Jz)
        push!(all_pz,pz)
    end
    indices = sortperm(all_pz);
    all_pz_sorted = all_pz[indices];
    aux = (scatter([1:tc],[all_pz_sorted],color=cor2,legend=false))
    display(aux)
    if sv
        nf = nfig; savefig("/home/ademir/Ademir/artigos/Mael/teste_quadprog$nf.pdf");
    end
    return all_pz
end
