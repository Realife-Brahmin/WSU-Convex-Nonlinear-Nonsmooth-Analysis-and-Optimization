using Parameters

include("../helperFunctions.jl")
include("../types.jl")

function updateTRModelParams(SR1params::Dict,
    # k::Int,
    # xk::Vector{Float64},
    # fk::Float64,
    # gk::Vector{Float64}
    solState::SolStateType
    )

    @unpack k, xk, fk, gk = solState
    n = length(xk)
    if k == 1
        Bk = fk*I(n)
    else
        @unpack xkm1, gkm1, Bkm1 = SR1params

        @checkForNaN sk = xk - xkm1
        @checkForNaN yk = gk - gkm1
        @checkForNaN wk = yk - Bkm1*sk

        ρinv = wk'*sk
        if ρinv == 0 || isnan(ρinv)
            @error "Bad ρ"
        else
            ρ = 1.0/ρinv
            Bk = Bkm1 + ρ*wk*wk'
        end
    end

    xkm1, gkm1, Bkm1 = xk, gk, Bk
    @pack! SR1params = xkm1, gkm1, Bkm1
    return SR1params
end