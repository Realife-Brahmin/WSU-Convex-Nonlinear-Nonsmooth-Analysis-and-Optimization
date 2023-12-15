using LinearAlgebra
using Parameters
using Test

include("helperFunctions.jl")
include("types.jl")


function getCandidateStep(SR1params::Dict,
    TRparams::Dict;
    itrMax::Int = 50)

    @unpack Delta = TRparams
    gk = SR1params[:gkm1] # it is now gk from updateModelParams()
    eps_k = min(1/2, sqrt(norm(gk)))*norm(gk)

    Bk = SR1params[:Bkm1] # Bkm1 was updated with Bk after the last updateTRModelParams

    j = 1
    keepFindingCandidate = true

    rj = gk
    dj = -rj
    zj = myzeros(gk)

    while keepFindingCandidate
        @checkForNaN ρinv = dj'*Bk*dj
        if ρinv ≤ 0
            pj, alphaj = getBoundaryIntersection(zj, dj, Delta)
            keepFindingCandidate =  false
            return pj
        else
            ρ = 1.0/ρinv
            alphaj = ρ*rj'*rj
            zjp1 = zj + alphaj*dj
        end

        if norm(zjp1) ≥ Delta # We've stepped outside the TRegion, so let's backup along this step to the boundary
            pj, alphaj = getBoundaryIntersection(zj, dj, Delta)
            keepFindingCandidate = false
            return pj
        else
            zj = zjp1
            rjp1 = rj + alphaj*Bk*dj
        end

        if norm(rjp1) < eps_k
            pj = zj
            keepFindingCandidate = false
            return pj
        else
            ρinv2 = rj'*rj
            if ρinv2 == 0 || isnan(ρinv2)
                @error "Invalid ρinv2"
            else
                ρ2 = 1.0/ρinv2
                betaj = ρ2*rjp1'*rjp1
                dj = -rjp1 + betaj*dj
            end
        end

        j += 1

        if j > itrMax
            @error "Failed to find TRegion step despite $(itrMax) iterations."
            keepFindingCandidate = false
        end

    end

    @error "floc"
end

"""
    getBoundaryIntersection(v1::Vector{Float64}, v2::Vector{Float64}, Delta::Float64)

Calculate the intersection point of the line defined by `v1 + alpha * v2` with the boundary of a trust region of radius `Delta`.

# Arguments
- `v1::Vector{Float64}`: Starting point of the vector.
- `v2::Vector{Float64}`: Direction vector.
- `Delta::Float64`: Radius of the trust region.

# Returns
- `(p, alpha)`: A tuple where `p` is the intersection point on the boundary of the trust region, and `alpha` is the scalar multiplier for `v2`.

# Errors
- Throws an error if `v1` and `v2` have mismatched lengths.
- Throws an error if `ρinv` (inverse of the dot product of `v2` with itself) is zero or NaN, indicating a problematic calculation.

# Example
```julia
v1 = [1.0, 2.0]
v2 = [3.0, 4.0]
Delta = 5.0
intersection = getBoundaryIntersection(v1, v2, Delta)
```
"""
function getBoundaryIntersection(v1::Vector{Float64},
    v2::Vector{Float64},
    Delta::Float64)

    if length(v1) != length(v2)
        @error "Cannot perform interpolation if vectors have mismatched lengths"
    else
        ρinv = v2'*v2
        if isnan(ρinv) || ρinv == 0
            @error "problematic ρinv"
        else
            ρ = 1.0/ρinv
            alpha = -ρ*v1'*v2 + sqrt( (ρ*v1'*v2)^2 - ρ*(v1'*v1 - Delta^2))
            p = v1 + alpha*v2
            return (p=p, alpha=alpha)
        end
    end

    @error "floc"
end