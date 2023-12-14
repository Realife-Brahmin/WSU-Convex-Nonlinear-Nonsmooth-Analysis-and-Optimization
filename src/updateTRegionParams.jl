using LinearAlgebra
using Parameters
using Test

include("helperFunctions.jl")
include("types.jl")

function updateTRegionParams(TRparams, ρk, pk, y)
    @unpack Delta, Delta_min, Delta_max, etta_1, etta_2, etta_3, delta_1, delta_2, updateRadius, accept = TRparams

    if ρk > etta_3 && norm(ρk) >= 0.99*Delta
        updateRadius = "Increase"
        Delta = min.(delta2*Delta, Delta_max)
        accept = true
        y += pk
        return y, TRparams
    elseif ρk > etta_2
        updateRadius = "Same"
        accept = true
        y += pk
    elseif ρk > etta_1
        updateRadius = "Decrease"
        Delta = max.(delta_1*Delta, Delta_min)
        accept = true
        y += pk
    else
        updateRadius = "Decrease"
        Delta = max.(delta_1*Delta, Delta_min)
        accept = false
    end

    @pack! TRparams = Delta, accept, updateRadius
    return y, TRparams

end