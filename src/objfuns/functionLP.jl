include("objective.jl")
include("../solveLP.jl")

using HiGHS
using JuMP
using Parameters

function LPObjectiveFunction(x::Vector{Float64}, 
    pDict;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)
    
    n = length(x)
    @unpack c = pDict
    f = transpose(x)*c

    if !getGradientToo
        return f
    elseif getGradientToo
        g = c
        return f, g
    else
        @error "floc"
    end
    
end
