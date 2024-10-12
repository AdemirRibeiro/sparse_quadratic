function quadprog_max2(n,s,m,η,c,V,λ,A,b,Jz)
   
   model = JuMP.Model(Ipopt.Optimizer)
   set_silent(model)

   @variable(model, α[1:n])
   @variable(model, β[1:m] ≥ 0)
   γ = @expression(model, c + V*(λ.*α) + A'*β)   
   @objective(model, Min, dot(b,β) + .25*dot(α,λ.*α) + .25*η*dot(γ[Jz],γ[Jz]))

   # print(model)
   optimize!(model);
   α, β = JuMP.value.(α), JuMP.value.(β)
   γ = c + V*(λ.*α) + A'*β;
   pz = -(dot(b,β) + .25*dot(α,λ.*α) + .25*η*dot(γ[Jz],γ[Jz]));
   return pz
end
