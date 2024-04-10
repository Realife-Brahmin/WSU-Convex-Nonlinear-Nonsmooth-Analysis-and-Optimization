using Base.Threads
using DataFrames
using LinearAlgebra
using Parameters

include("../helperFunctions.jl");
include("boxConstraintPenalty.jl");
include("equalityConstrainedQP.jl");
include("../../alg.jl"); # include alg.jl from parent directory
include("../types.jl")

function generate_pr(functionName::Function,
    x0::Vector{Float64};
    problemType::String="Unconstrained",
    method::String="QuasiNewton",
    data=Matrix{Float64}(undef, 0, 0),
    params=Dict(),
    objectiveString::String="undefinedFunction",
    verbose::Bool=true,
)

    p = Dict(:params => params, :data => data)

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