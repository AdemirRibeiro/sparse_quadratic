function quadprog_max(n,s,m,η,c,V,λ,A,b,Jz)
   
   sqrtλ = sqrt.(λ);
   model = JuMP.Model(Ipopt.Optimizer)
   set_silent(model)

   @variable(model, α[1:n])
   @variable(model, β[1:m] ≥ 0)
   γ = @expression(model, c + V*(sqrtλ.*α) + A'*β)
   @objective(model, Min, dot(b,β) + .25*dot(α,α) + .25*η*dot(γ[Jz],γ[Jz]))

   # print(model)
   optimize!(model);
   α, β = JuMP.value.(α), JuMP.value.(β)
   γ = c + V*(sqrtλ.*α) + A'*β;
   pz = -(dot(b,β) + .25*dot(α,α) + .25*η*dot(γ[Jz],γ[Jz]));
   return pz
end
