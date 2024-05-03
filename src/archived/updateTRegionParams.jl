using LinearAlgebra
using Parameters
using Test

include("../helperFunctions.jl")
include("../types.jl")

function updateTRegionParams(
    TRparams, 
    ρk, 
    pk, 
    y;
    verbose::Bool = true)

    @unpack Delta, Delta_min, Delta_max, etta_1, etta_2, etta_3, delta_1, delta_2, updateRadius, accept = TRparams

    if ρk > etta_3 && norm(pk) >= 0.99*Delta
        updateRadius = "Increase"
        updateRadiusString = "increased"
        Delta = min.(delta_2*Delta, Delta_max)
        accept = true
        y += pk
        return y, TRparams
    elseif ρk > etta_2
        updateRadius = "Same"
        updateRadiusString = "the same"
        accept = true
        y += pk
    elseif ρk > etta_1
        updateRadius = "Decrease"
        updateRadiusString = "decreased"
        Delta = max.(delta_1*Delta, Delta_min)
        accept = true
        y += pk
    else
        updateRadius = "Decrease"
        updateRadiusString = "decreased"
        Delta = max.(delta_1*Delta, Delta_min)
        accept = false
    end

    accept == true ? acceptString = "accepted" : acceptString = "rejected"
    
    myprintln(verbose, "Candidate Step: y = $(y) has been $(acceptString) and the trust region has been set to be $(updateRadiusString)")
    
    @pack! TRparams = Delta, accept, updateRadius
    return y, TRparams

end