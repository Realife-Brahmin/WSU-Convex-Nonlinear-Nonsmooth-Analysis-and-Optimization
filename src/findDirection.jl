include("helperFunctions.jl")

# For ConjugateGradientDescent update
mutable struct CGargsType
    k::Int
    xk::Vector{Float64}
    xkp1::Vector{Float64}
    gk::Vector{Float64}
    gkp1::Vector{Float64}
end

function constructorGCargs(
    pr::NamedTuple)::CGargsType
    k = 1
    xk = pr.x0
    xkp1 = similar(xk)
    gk = similar(xk)
    gkp1 = similar(xk)
    CGargs = CGargsType(k, xk, xkp1, gk, gkp1)
    return CGargs
end

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
    k = 1
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
    QNargs::QNargsType=constructorQNargs(pr),
    CGargs::CGargsType=constructorGCargs(pr),
    verbose::Bool=false)

    method = pr.alg.method
    n = length(∇fk)
    if method == "GradientDescent"
        Bₖ = I(n)
    elseif method == "ConjugateGradientDescent"
        @show k = CGargs.k
        @show CGargs
        if k == 1
            ∇f0 = ∇fk
            pkp1 = -∇f0 
        else
            ∇fkp1 = CGargs.gkp1
            pk = CGargs.pk
            βkp1 = max(0, ∇fkp1'*(∇fkp1-∇fk)/(∇fkp1'*∇fk))
            pkp1 = -∇fkp1 + βkp1*pk
        end
        pₖ = pkp1
        GCargs.pk = pₖ
        CGargs.xk = xkp1
        CGargs.gk = gkp1
        return pₖ, CGargs

    elseif method == "QuasiNewton"
        @unpack k, xk, xkp1, fk, gk, gkp1, Hk = QNargs
        @show k
        if k == 1
            H0 = QNargs.fk * I(n)
            QNargs.Hk = H0
            QNargs.xk = xkp1
            QNargs.gk = gkp1
            Bₖ = H0
        else
            # @checkForNaN(gk)
            # @show gk
            # @checkForNaN(gkp1)
            # @show gkp1
            # @checkForNaN(xk)
            # @checkForNaN(xkp1)
            sk = xkp1 - xk
            yk = gkp1 - gk
            ρkinv = yk'*sk
            if ρkinv == 0
                @warn "So, Hkp1 is actualy the same as Hk?, Maybe stop this?"
            end
            ρk = 1.0/ρkinv
            if ρk < 0
                @warn "QuasiNewton step problematic! y'*s < 0!"
                myprintln(true, "Making a different H from current value of f.")
                Hkp1 = fk*I(n)
            elseif isnan(ρk)
                @error "NaN!"
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