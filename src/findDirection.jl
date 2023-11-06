include("helperFunctions.jl")
include("types.jl")

function findDirection(
    pr::NamedTuple, gk::Vector{Float64};
    CGState::CGStateType=CGStateType(),
    QNState::QNStateType=QNStateType(),
    verbose::Bool=false)

    n = length(gk)
    method = pr.alg.method
    
    if method == "GradientDescent"
        Bₖ = I(n)
        pₖ = -Bₖ*gk
        return pₖ

    elseif method == "ConjugateGradientDescent"
        @unpack k, gk, pk, justRestarted = CGState

        @show k
        if k == 1
            pk = -gk 
            betak = 0.0
        else
            @unpack gkm1, pkm1 = CGState
            diff = gk'*(gk-gkm1)
            mag = gkm1'*gkm1
            @show betak = max(0, diff/mag)
            if betak == 0
                justRestarted = true
                myprintln(true, "Restarted ConjugateGradientDescent.")
            end
            @checkForNaN pk = -gk + betak*pkm1
        end

        gkm1 = gk
        pkm1 = pk
        betakm1 = betak
        @pack! CGState = gkm1, pkm1, betakm1
        pₖ = pk
        return pₖ, CGState

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
        pₖ = -Bₖ*gk

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

    @error "floc"

end