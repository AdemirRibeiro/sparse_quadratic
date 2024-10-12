function quadprog_max3(n,s,m,η,c,V,λ,A,b,Jz)
   Jlb = findall(λ .> 0);
   λJ = λ[Jlb]; 
   Ilb = setdiff(collect(1:n), Jlb);
   zeros_Ilb = 0*Ilb; 
   
   model = JuMP.Model(Ipopt.Optimizer)
   set_silent(model)

   @variable(model, α[1:n])
   @variable(model, β[1:m] ≥ 0)
   γ = @expression(model, c + V*α + A'*β)  
   αJ = @expression(model, α[Jlb])  
   @objective(model, Min, dot(b,β) + .25*dot(αJ,(αJ./λJ)) + .25*η*dot(γ[Jz],γ[Jz]))
   @constraint(model, α[Ilb] .== zeros_Ilb)

   # print(model)
   optimize!(model);
   α, β = JuMP.value.(α), JuMP.value.(β)
   γ = c + V*α + A'*β;   
   αJ = α[Jlb];  
   pz = -(dot(b,β) + .25*dot(αJ,(αJ./λJ)) + .25*η*dot(γ[Jz],γ[Jz]));
   return pz
end


