
using DataFrames
# using GR
using LaTeXStrings
using Parameters
using Plots
using Printf

# include("utilities.jl")
include("objfuns/objective.jl")

function showresults(res::NamedTuple;
    log::Bool=true,
    log_path::String="./logging/")

    @unpack pr = res
    @unpack method = pr[:alg]

    if method == "NelderMead"
        @unpack converged, statusMessage, xvals, fvals, fevals, cause = res
        gevals = 0
        result_txt = log_path * "results_" * string(pr.objective) * "_" * pr.alg[:method] * "_" * string(pr.alg[:maxiter]) * ".txt"
        lsMsg = ""

        tolMsg = "Simplex Diameter Tolerance = $(res.pr[:alg][:DeltaTol])"

    else
        @unpack converged, statusMessage, fvals, αvals, backtrackVals, xvals, M, fevals, gevals, cause = res
        result_txt = log_path * "results_" * string(pr.objective) * "_" * pr.alg[:method] * "_" * pr.alg[:linesearch] * "_" * string(pr.alg[:maxiter]) * ".txt"
        lsMsg = "Linesearch method used: " * pr.alg[:linesearch]

        tolMsg = "dftol = $(res.pr.alg[:dftol])"

    end

    params = pr.p[:params]

    if isfile(result_txt)
        rm(result_txt) # remove results log file if already present
    end
    log, log_path = log, result_txt  # or whatever default values you have

    function myprintln1(v, message)
        myprintln(v, message, log=log, log_path=result_txt)
    end

    
    nnztol = 1e-8 # variable having a value below this will NOT be counted in the list of non-zero variables. For printing purposes only. 

    v = true
    myprintln1(v, "****************************")
    myprintln1(true, "Solver run has concluded.")
    myprintln1(true, "Function solved for: $(pr.objective)()")
    myprintln1(true, "Cause(s) for stopping:")
    causeForStopping = res.cause
    for msg ∈ causeForStopping
        myprintln1(true, msg)
    end
    if converged == true
        myprintln1(true, statusMessage)
    elseif converged == false
        myprintln1(true, statusMessage)
    else
        @error "bad condition"
    end

    methodMsg = "Method Used: " * pr.alg[:method]
    myprintln1(v, "****************************")
    myprintln1(v, methodMsg)
    myprintln1(v, lsMsg)

    
    n = size(xvals, 1)
    itr = size(xvals, 2)
    x☆ = xvals[:, itr]

    myprintln1(v, "****************************")
    
    if converged == true
        fvalPrefixMinimum = "Optimal function value = "
        xvalMessage = "Optimal Nonzero Variables:"
    elseif converged == false
        fvalPrefixMinimum = "Best function value = "
        xvalMessage = "Best Known Nonzero Variables:"
    else
        @error "bad condition"
    end

    fvalMessageMinimum = fvalPrefixMinimum*"$(fvals[itr])"

    myprintln1(v, fvalMessageMinimum)

    myprintln1(v, "***************************")
    myprintln1(v, xvalMessage)
    if n <= 15
        # If the vector is short enough, print all elements
        for k in 1:n
            myprintln1(v, "x☆[$k] = $(x☆[k])")
        end
    else
        # If the vector is long, print an even distribution of elements
        myprintln1(v, "Only upto 15 values will be printed.")

        mid_indices = round.(Int, LinRange(6, n-5, 7)[2:end-1])
        indices_to_print = unique([1:5; mid_indices; n-4:n])
        for k in indices_to_print
            myprintln1(v, "x☆[$k] = $(x☆[k])")
        end
    end

    myprintln1(v, "***************************")
    myprintln1(v, "Number of fevals = $(fevals)")
    myprintln1(v, "Number of gevals = $(gevals)")
    myprintln1(v, "Number of evals = $(fevals+gevals*n)")

    myprintln1(v, "*"^50)
    myprintln1(v, "Tolerance critera used:")
    myprintln1(v, tolMsg)
    
end

