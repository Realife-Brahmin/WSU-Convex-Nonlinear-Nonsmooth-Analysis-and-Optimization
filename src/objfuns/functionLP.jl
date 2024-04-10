include("objective.jl")

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

pLP = Dict(:c => c, :Ae => Ae, :be => be, :Ai => Ai, :bi => bi)

objective = LPObjectiveFunction
objectiveOriginal = LPObjectiveFunction
objectiveString = string(objectiveOriginal)
params = pLP

pr = generate_pr(objective, w0, params=params, problemType="LP"; objectiveString=objectiveString)

# f0 = equalityConstrainedQP(w0, pECQP)

# solState = SolStatePGCGType(w0, G, c, A, fk=f0)
