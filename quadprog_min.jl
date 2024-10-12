function quadprog_min(n,s,m,η,c,V,λ,A,b,Jz)

   sqrtλ = sqrt.(λ);
   Iz=setdiff(collect(1:n), Jz);
   zeros_Iz=0*Iz;

   model = JuMP.Model(Ipopt.Optimizer)
   set_silent(model)

   @variable(model, x[1:n])
   @variable(model, y[1:n])
   # lby2 = @expression(model, λ.*(y.^2))
   # @objective(model, Min, dot(c,x) + sum(lby2) + (1/η)*dot(x,x))
   @objective(model, Min, dot(c,x) + dot(λ,(y.^2)) + (1/η)*dot(x,x))
   @constraint(model, A*x .≤ b)
   # @constraint(model, [j=1:n], y[j] - dot(V[:,j],x) == 0)
   @constraint(model, [j=1:n], sqrtλ[j]*y[j] - sqrtλ[j]*dot(V[:,j],x) == 0)
   @constraint(model, x[Iz] .== zeros_Iz)

   # print(model)
   optimize!(model)
   x, y  = JuMP.value.(x), JuMP.value.(y)
   lby2 = λ.*(y.^2)
   pz = dot(c,x) + sum(lby2) + (1/η)*dot(x,x);
   return pz
end
