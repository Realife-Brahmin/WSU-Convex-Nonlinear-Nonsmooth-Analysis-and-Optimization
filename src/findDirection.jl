# For QasiNewton BFGS update
mutable struct QNargsType
    k::Int
    xk::Vector{Float64}
    xkp1::Vector{Float64}
    fk::Float64
    gk::Vector{Float64}
    gkp1::Vector{Float64}
    Hk::Matrix{Float64}
end


function constructorQNargs(
    pr::NamedTuple; 
    fk=pr.objective(pr.x0, pr.p, getGradientToo=false))::QNargsType
    k = 0
    xk = pr.x0
    n = length(xk)
    xkp1 = similar(xk)
    gk = similar(xk)
    gkp1 = similar(xk)
    Hk = fk*I(n)
    QNargs = QNargsType(k, xk, xkp1, fk, gk, gkp1, Hk)
    return QNargs
end

function findDirection(
    pr::NamedTuple, ∇fk::Vector{Float64};
    QNargs::QNargsType=constructorQNArgs(pr),
    verbose::Bool=false)

    method = pr.alg.method
    n = length(∇fk)
    if method == "GradientDescent"
        Bₖ = I(n)
    elseif method == "ConjugateGradientDescent"
        @error "Currently not formulated for this method"
    elseif method == "QuasiNewton"
        k = QNargs.k
        if k == 1
            H0 = QNargs.fk * I(n)
            Bₖ = H0
        else
            @unpack xk, xkp1, fk, gk, gkp1, Hk = QNargs
            sk = xkp1 - xk
            yk = gkp1 - gk
            ρk = 1.0/(yk'*sk)
            if ρk < 0
                @warn "QuasiNewton step problematic! y'*s < 0!"
                myprintln(true, "Making a different H from current value of f.")
                Hkp1 = fk*I(n)
            else
                Hkp1 = (I(n) - ρk*sk*yk')*Hk'*(I-ρk*yk*sk') + ρk*sk*sk'
            end

            QNargs.Hk = Hkp1
            QNargs.xk = xkp1
            QNargs.gk = gkp1
            Bₖ = Hkp1
        end
        pₖ = -Bₖ*∇fk
        return pₖ, QNargs
    else
        @error "Currently not formulated for this method"
    end
    pₖ = -Bₖ*∇fk
    return pₖ
end