function prog_maxmin2(n,s,m,η,c,V,λ,A,b)

   α = zeros(n);
   β = zeros(m);
   γ = c + V*(λ.*α) + A'*β;
   indices = sortperm(abs.(γ),rev = true);
   Jz = indices[1:s];
   ξstar = dot(b,β) + .25*dot(α,λ.*α) + .25*η*dot(γ[Jz],γ[Jz]);
   all_ξstar = Vector{Any}();
   push!(all_ξstar,ξstar)
   k = 0;
   k_max = 5000;
   r = 1/2;
   tol = 1e-5;
   for j = 0:k_max
      k = k + 1;
      # Computing diag(z)*γ, denoted by zγ
      Iz = setdiff(collect(1:n), Jz);
      zeros_Iz = 0*Iz;
      aux = γ;
      aux[Iz] = zeros_Iz;
      zγ = aux;
      g_α = .5*λ.*α + .5*η*diagm(λ)*(V')*zγ;
      g_β = b + .5*η*A*zγ;
      ng_α = norm(g_α);
      ng_β = norm(g_β);
      ng = sqrt(ng_α^2 + ng_β^2);
      tk = 1/(k^r);
      α = α - tk*g_α/ng;
      β = max.(zeros(m),β - tk*g_β/ng);
      γ = c + V*(λ.*α) + A'*β;
      indices = sortperm(abs.(γ),rev=true);
      Jz = indices[1:s];
      ξstar_old = ξstar;
      ξ = dot(b,β) + .25*dot(α,λ.*α) + .25*η*dot(γ[Jz],γ[Jz]);
      ξstar = min(ξstar_old, ξ);
      push!(all_ξstar,ξstar)
      diffξstar = abs(ξstar-ξstar_old)/(abs(ξstar) + 0.00001);
      # if diffξstar < tol
      #    break
      # end
   end
   d_opt = -ξstar;
   return d_opt#, all_ξstar[1:8], all_ξstar[end-7:end]
end

