function data(num_prob)
    all_n = rand(MersenneTwister(15),6:15,num_prob);
    all_m = rand(MersenneTwister(45),2:6,num_prob);
    η1 = rand(MersenneTwister(56),1:3,num_prob);
    η2 = rand(MersenneTwister(65),1:10,num_prob);
    all_η = η1./η2;
    all_s = Vector{Int64}();
    all_c = Vector{Any}();
    all_V = Vector{Any}();
    all_λ = Vector{Any}();
    all_A = Vector{Any}();
    all_b = Vector{Any}();
    s1 = Int.(ceil.(all_n/4));
    s2 = Int.(ceil.(all_n/2));
    aux_seed = rand(1:1000,6);
    aux_seed = [580, 590, 613, 241, 552, 806]
    # aux_seed = [919, 299, 480, 761, 823, 432];
    # aux_seed = [259  717  805  952  371  478];
    all_ns1 = rand(MersenneTwister(aux_seed[1]),1:1000,num_prob);
    all_ns2 = rand(MersenneTwister(aux_seed[2]),10:1100,num_prob);
    all_ns3 = rand(MersenneTwister(aux_seed[3]),20:2000,num_prob);
    all_ns4 = rand(MersenneTwister(aux_seed[4]),142:1700,num_prob);
    all_ns5 = rand(MersenneTwister(aux_seed[5]),131:10000,num_prob);
    all_ns6 = rand(MersenneTwister(aux_seed[6]),211:1080,num_prob);
    for np = 1:num_prob
        ns1 = all_ns1[np]; seed1 = MersenneTwister(ns1);
        ns2 = all_ns2[np]; seed2 = MersenneTwister(ns2);
        ns3 = all_ns3[np]; seed3 = MersenneTwister(ns3);
        ns4 = all_ns4[np]; seed4 = MersenneTwister(ns4);
        ns5 = all_ns5[np]; seed5 = MersenneTwister(ns5);
        ns6 = all_ns6[np]; seed6 = MersenneTwister(ns6);
        s = rand(seed1,s1[np]:s2[np]);
        push!(all_s,s)
        n = all_n[np];
        m = all_m[np];
        c = 5*(rand(seed2,n).-.5);
        push!(all_c,c)
        aux = 5*(rand(seed3,n,n).-.5);
        P, R = qr(aux);
        V = Matrix(P);
        push!(all_V,V)
        λ = 5*(rand(seed4,n));
        # Set some lambdas to zero or very small values
        qt_lb_small = Int(ceil(n/6));
        λ[1:qt_lb_small] = zeros(qt_lb_small);
        λ[end-qt_lb_small+1:end] = (1e-10)*ones(qt_lb_small);
        push!(all_λ,λ)
        A = 5*(rand(seed5,m,n).-.5);
        push!(all_A,A)
        b = 5*rand(seed6,m);
        # b = 15*(rand(seed6,m).-.5);
        push!(all_b,b)
    end
    return aux_seed, all_n, all_s, all_m, all_η, all_c, all_V, all_λ, all_A, all_b
end
