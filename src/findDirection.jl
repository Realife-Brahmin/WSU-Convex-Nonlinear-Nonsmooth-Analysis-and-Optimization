function constructorQNargs(;)
    k = 0
    xk = pr.x0
    xkp1 = similar(xk)
    gk = similar(xk)
    gkp1 = similar(xk)
    Hk = pr.objective(pr.x0, pr.p, getGradientToo=false)
    QNargs = (k=k, xk=xk, xkp1=xkp1, gk=gk, gkp1=gkp1, Hk=Hk)
    return QNargs
end

function findDirection(
    pr::NamedTuple, ∇fk::Vector{Float64};
    QNargs::NamedTuple=constructorQNArgs(),
    verbose::Bool=false)::Vector{Float64}

    method = pr.alg.method
    n = length(∇fk)
    if method == "GradientDescent"
        Bₖ = I(n)
        # # pₖ = -Bₖ*∇fk
        # pₖ = -∇fk
    elseif method == "ConjugateGradientDescent"
        @error "Currently not formulated for this method"
    elseif method == "QuasiNewton"
        k = QNargs.k
        if k == 0
            H0 = QNargs.Hk
            Bₖ = H0
        else
            @unpack xk, xkp1, ∇fk, ∇fkp1, Hk = QNargs
            sk = xkp1 - xk
            yk = ∇fkp1 - ∇fk
            ρk = 1.0/(yk'*sk)
            if ρk < 0
                @warn "QuasiNewton step problematic! y'*s < 0!"
            end
            Hkp1 = (I(n) - ρk*sk*yk')*Hk'*(I-ρk*yk*sk') + ρk*sk*sk'
            Bₖ = Hkp1
        end
    else
        @error "Currently not formulated for this method"
    end
    pₖ = -Bₖ*∇fk
    return pₖ
end

function constructorQNArgs()
end