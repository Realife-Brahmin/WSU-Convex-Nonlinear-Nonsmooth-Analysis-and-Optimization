include("helperFunctions.jl")
include("types.jl")

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

function findDirection(
    pr::NamedTuple, ∇fk::Vector{Float64};
    CGargs::CGargsType=constructorCGargs(pr),
    CGState::CGStateType=CGStateType(),
    QNState::QNStateType=QNStateType(),
    verbose::Bool=false)

    method = pr.alg.method
    n = length(∇fk)
    

    gk = ∇fk

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
        @unpack k, xkm1, xk, fkm1, fk, gkm1, gk, Hkm1, Hk = QNState

        @show k
        if k == 1
            H0 = fk * I(n)
            Hk = H0
        else
            sk = xk - xkm1
            yk = gk - gkm1
            ρkinv = yk'*sk
            if ρkinv == 0
                @warn "So, Hk is actually the same as Hkm1?, Maybe stop this?"
            end
            ρk = 1.0/ρkinv
            if ρk < 0
                @warn "QuasiNewton step problematic! y'*s < 0!"
                myprintln(true, "Making a different H from current value of f.")
                Hk = fk*I(n)
            elseif isnan(ρk)
                @error "NaN!"
            else
                Hk = (I(n) - ρk*sk*yk')*Hkm1'*(I-ρk*yk*sk') + ρk*sk*sk'
            end
        end

        Bₖ = Hk
        pₖ = -Bₖ*∇fk

        Hkm1 = Hk
        xkm1 = xk
        fkm1 = fk
        gkm1 = gk
        pkm1 = pₖ
        
        @pack! QNState = xkm1, fkm1, gkm1, pkm1, Hkm1
        return pₖ, QNState

    else
        @error "Currently not formulated for this method"
    end

    pₖ = -Bₖ*∇fk
    return pₖ
end