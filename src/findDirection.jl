include("helperFunctions.jl")

# For ConjugateGradientDescent update
mutable struct CGargsType
    k::Int
    gk::Vector{Float64}
    # gkp1::Vector{Float64}
    pk::Vector{Float64}
    justRestarted::Bool
end

function constructorCGargs(
    pr::NamedTuple)::CGargsType
    k = 1
    xk = pr.x0
    gk = myzeros(xk)
    # gkp1 = myzeros(xk)
    pk = myzeros(xk)
    justRestarted = false
    CGargs = CGargsType(k, gk, pk, justRestarted)
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
    xkp1 = myzeros(xk)
    gk = myzeros(xk)
    gkp1 = myzeros(xk)
    Hk = fk*I(n)
    QNargs = QNargsType(k, xk, xkp1, fk, gk, gkp1, Hk)
    return QNargs
end

function findDirection(
    pr::NamedTuple, ∇fk::Vector{Float64};
    QNargs::QNargsType=constructorQNargs(pr),
    CGargs::CGargsType=constructorCGargs(pr),
    CGState::CGStateType=CGStateType(),
    verbose::Bool=false)

    method = pr.alg.method
    n = length(∇fk)
    gkp1 = ∇fk
    if method == "GradientDescent"
        Bₖ = I(n)
    elseif method == "ConjugateGradientDescent"
        # @unpack k, gk, pk, justRestarted = CGargs
        @unpack k, fk, gk, gkmag, pk = CGState

        @show k
        if k == 1
            pkp1 = -gkp1 
        else
            diff = gkp1'*(gkp1-gk)
            mag = gk'*gk
            βkp1 = max(0, diff/mag)
            if βkp1 == 0
                justRestarted = true
                myprintln(true, "Restarted ConjugateGradientDescent.")
            else
                justRestarted = false
            end
            @checkForNaN pkp1 = -gkp1 + βkp1*pk
        end

        CGargs.pk = pkp1
        CGargs.gk = gkp1
        CGargs.justRestarted = justRestarted

        pₖ = pkp1
        return pₖ, CGargs

    elseif method == "QuasiNewton"
        @unpack k, xk, xkp1, fk, gk, gkp1, Hk = QNargs
        @show k
        if k == 1
            H0 = QNargs.fk * I(n)
            Hkp1 = H0
        else
            sk = xkp1 - xk
            yk = gkp1 - gk
            ρkinv = yk'*sk
            if ρkinv == 0
                @warn "So, Hkp1 is actually the same as Hk?, Maybe stop this?"
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

            # Bₖ = Hkp1
        end

        Bₖ = Hkp1
        QNargs.Hk = Hkp1
        QNargs.xk = xkp1
        QNargs.gk = gkp1

        pₖ = -Bₖ*∇fk
        return pₖ, QNargs
    else
        @error "Currently not formulated for this method"
    end
    pₖ = -Bₖ*∇fk
    return pₖ
end