using LinearAlgebra, JuMP, Ipopt, Random, Statistics, DataFrames, CSV, Distributions


function CreateData_Unfair_SVM(np, seed)
    #Gera os betas, os b's e os X randomicos e aí acha as labels...
    kkkk1 = np
    #β = [-5;4;8;5;6]./10
    β = [-20;4;8;5;20]./10

    #b =  rand(MersenneTwister(seed), Normal(0, 3.0), nc)
    sexH = rand(MersenneTwister(seed), kkkk1, 1) .< 0.5
    #race = rand(MersenneTwister(seed+1), kkkk1, 1) .< 0.5
    #marital = rand(MersenneTwister(seed+2), kkkk1, 1) .< 0.5

    mean = [0.0, 0.0, 0.0]
    C = [1.0 0 0; 0 1.0 0; 0 0 1.0]
    d = MvNormal(mean, C)
    #X = [ones(kkkk1) rand(MersenneTwister(seed), d, kkkk1)' race marital sexH]
    X = [ones(kkkk1) rand(MersenneTwister(seed), d, kkkk1)' sexH]


    function predicaoSVM(β,X)
        

        m12(x,β) = dot(β,x)
        ly = size(X)[1]
        verificar = zeros(ly)
        y = zeros(ly)

        for i = 1 : ly
            a = m12(X[i,:],β)
            verificar[i] = exp(a)/(1+exp(a))
            d = [Binomial(1, verificar[i])]
            y[i] = rand(MersenneTwister(42),d[1], 1)[1]

            if y[i] == 0
                y[i] = -1
            elseif y[i] == 1
                y[i] = 1
            end

        end
        return y
    end
    y = predicaoSVM(β,X)

    target = zeros(kkkk1)
    for i = 1 : kkkk1
        if y[i] == 1
            target[i] = 1
        else
            target[i] = -1
        end
    end

    sex = sexH
    CDT = [X[:,2:end] target]

    CDT = DataFrame(CDT, :auto)
    CDT = rename(CDT, :x1 => :F1,:x2 => :F2,:x3 => :F3,:x4 => :Sex, :x5 => :target)#, :x5 => :Race, :x6 => :Marital, :x7 => :target)

    CDT = shuffle(MersenneTwister(42), CDT)
    sizeXTrain = 0.01
    df1 = CDT[1:Int64(round(sizeXTrain*size(CDT,1))),:]
    df2 = CDT[Int64(round(sizeXTrain*size(CDT,1)))+1:end,:]

    X_train = df1[:,1:end-1]
    Y_train = df1[:,end]
    New_data = df2[:,1:end-1]
    Y_New_data = df2[:,end]

    dif = rand(MersenneTwister(42), 1:length(Y_train), Int64(round(0.2*size(Y_train,1))))
    for i in dif
        if Y_train[i] == 1
            Y_train[i] = -1
        else
            Y_train[i] = 1
        end
    end

    dif = rand(MersenneTwister(42), 1:length(Y_New_data), Int64(round(0.2*size(Y_New_data,1))))
    for i in dif
        if Y_New_data[i] == 1
            Y_New_data[i] = -1
        else
            Y_New_data[i] = 1
        end
    end

    return X_train,Y_train,New_data,Y_New_data
end

X_train, Y_train, X_test, Y_test = CreateData_Unfair_SVM(10000, 42)



#A = rand(5,4)
#b = rand(5)

m,n = size(X_train)

model = JuMP.Model(Ipopt.Optimizer)

@variable(model, β[1:n])
@variable(model, ξ[1:m] ≥ 0)

@objective(model, Min, dot(β, β)/2 + sum(ξ))
@constraint(model, [i=1:m], (dot(β, X_train[i,:])) * Y_train[i] ≥ 1 - ξ[i])
#@constraint(model,  A*β .≤ b)

print(model)
optimize!(model)
β, ξ  = JuMP.value.(β), JuMP.value.(ξ)
