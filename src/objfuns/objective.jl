using Base.Threads
using DataFrames
using LinearAlgebra
using Parameters

include("../helperFunctions.jl");
# include("../activeSetQP")
# include("boxConstraintPenalty.jl");
# include("../equalityConstrainedQP.jl");
include("../../alg.jl"); # include alg.jl from parent directory
include("../types.jl")


function generate_pr(functionName::Function,
    x0::Vector;
    problemType::String="Unconstrained",
    method::String="QuasiNewton",
    params=Dict(),
    objectiveString::String="undefinedFunction",
    verbose::Bool=true,
)
    
    p = params
    if objectiveString == "undefinedFunction"
        objectiveString = string(functionName)
    elseif objectiveString != string(functionName)
        myprintln(verbose, "A wrapper function $(string(functionName)) will be used to solve for the original problem of $(objectiveString)")
    end

    alg = create_algorithm_settings(problemType=problemType, method=method)

    pr = (p=p, x0=x0, objective=functionName, alg=alg, problemType=problemType, objectiveString=objectiveString)

    myprintln(verbose, "Problem (pr::NamedTuple) generated")

    return pr
end