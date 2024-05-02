
using DataFrames
using LaTeXStrings
using Parameters
using Plots
using Printf

include("objfuns/objective.jl")

"""
    showresults(res; log::Bool=true, log_path::String="./logging/")

Generate detailed logs and outputs for optimization results, handling both single and multiple run scenarios.

# Arguments
- `res`: The result or vector of results from optimization processes. Each result should contain detailed optimization data including the problem setup, final solution, and diagnostic information.

# Keyword Arguments
- `log::Bool=true`: Indicates whether to log the results to a file.
- `log_path::String="./logging/"`: The directory path where result logs should be saved. This path is used to construct the filename for saving the results.

# Notes
- This function is designed to provide a comprehensive log of the optimization process, making it invaluable for detailed analysis and troubleshooting.
- Ensure the `log_path` directory exists or is created before running this function to avoid errors in file handling.
- The function is capable of handling various optimization methods and adjusts its output based on the specific method and settings used in each optimization result.
- Depending on the method, additional specific information like tolerance levels and line search details are reported to provide full insight into the optimization behavior.

# Description
This function iterates through the optimization results, processing and logging detailed information about each run. It supports both individual and batch optimization results. For each result, it extracts and logs detailed information including the optimization method used, convergence status, final and best function values, solution vectors, and other diagnostic statistics.

# Usage
- For single optimization results, it processes the contained information and logs it accordingly.
- For vectors containing multiple optimization results, it recursively calls itself to handle and log each individual result.

# Examples
```julia
# Assuming `results` is a single optimization result or a vector of results
showresults(results, log=true, log_path="./logs/")

# This will process and log detailed information about the optimization process for each result in `results`.
```
"""
function showresults(res;
    log::Bool=true,
    log_path::String="./logging/")

    if typeof(res) <: Vector
        numRuns = length(res)
        for runNum ∈ 1:numRuns
            res1 = res[runNum]
            showresults(res1)
        end

        return

    else
        
        res1 = res
    end

    @unpack pr = res1
    @unpack method = pr[:alg]
    objString = pr.objectiveString

    if method == "ActiveSetQP"
        @unpack converged, statusMessage, xvals, fvals, fevals, cause = res1
        @unpack etol, itol = pr[:alg]
        gevals = 0

        result_txt = log_path * "results_" * objString * "_" * pr.alg[:method] * "_" * string(pr.alg[:maxiter]) * ".txt"

        lsMsg = ""

        tolMsg = "PGCG Tolerance = $(etol)\n"*"ASQP Tolerance = $(itol)"

    elseif method == "AugmentedLagrangian"
        @unpack converged, statusMessage, xvals, fvals, fevals, cause = res1
        @unpack dxtol = pr[:alg]
        gevals = 0

        result_txt = log_path * "results_" * objString * "_" * pr.alg[:method] * "_" * string(pr.alg[:maxiter]) * ".txt"

        lsMsg = ""

        tolMsg = "dx Tolerance = $(dxtol)"

    elseif method == "NelderMead"
        @unpack converged, statusMessage, xvals, fvals, fevals, cause = res1
        gevals = 0

        result_txt = log_path * "results_" * objString * "_" * pr.alg[:method] * "_" * string(pr.alg[:maxiter]) * ".txt"

        lsMsg = ""

        tolMsg = "Simplex Diameter Tolerance = $(res1.pr[:alg][:DeltaTol])"

    elseif method == "GeneticAlgorithm"

        @unpack converged, statusMessage, xvals, fvals, fevals, cause = res1
        @unpack dftol = pr[:alg]
        gevals = 0

        result_txt = log_path * "results_" * objString * "_" * pr.alg[:method] * "_" * string(pr.alg[:maxiter]) * ".txt"

        lsMsg = ""

        tolMsg = "Generation Fitness Stagnation Tolerance = $(Int(res1.pr[:alg][:fvalRepeatTol])) consecutive iterations at dftol = $(dftol)."
    
    elseif method == "ProjectedGradientCG"
        @unpack converged, statusMessage, xvals, fvals, fevals, cause = res1
        @unpack etol = pr[:alg]
        gevals = 0

        result_txt = log_path * "results_" * objString * "_" * pr.alg[:method] * "_" * string(pr.alg[:maxiter]) * ".txt"

        lsMsg = ""

        tolMsg = "PGCG Tolerance = $(etol)"

    else
        @unpack converged, statusMessage, fvals, αvals, backtrackVals, xvals, fevals, gevals, cause = res1
        result_txt = log_path * "results_" * objString * "_" * pr.alg[:method] * "_" * pr.alg[:linesearch] * "_" * string(pr.alg[:maxiter]) * ".txt"
        lsMsg = "Linesearch method used: " * pr.alg[:linesearch]

        tolMsg = "dftol = $(res1.pr.alg[:dftol])"

    end

    # params = pr.p[:params]
    params = pr.p

    if isfile(result_txt)
        rm(result_txt) # remove results log file if already present
    end
    log, log_path = log, result_txt  # or whatever default values you have

    function myprintln1(v, message; color=:normal)
        myprintln(v, message, log=log, log_path=result_txt, color=color)
    end

    
    nnztol = 1e-8 # variable having a value below this will NOT be counted in the list of non-zero variables. For printing purposes only. 

    v = true
    myprintln1(v, "****************************")
    myprintln1(true, "Solver run has concluded.", color=:light_cyan)
    myprintln1(true, "Function solved for: $(objString)()", color=:light_cyan)
    myprintln1(true, "Cause(s) for stopping:", color=:light_cyan)
    causeForStopping = res1.cause
    for msg ∈ causeForStopping
        myprintln1(true, msg, color=:normal)
    end
    if converged == true
        myprintln1(true, statusMessage, color=:green)
    elseif converged == false
        myprintln1(true, statusMessage, color=:red)
    else
        @error "bad condition"
    end

    methodMsg = "Method Used: " * pr.alg[:method]
    myprintln1(v, "****************************")
    myprintln1(v, methodMsg)
    myprintln1(v, lsMsg)

    
    n = size(xvals, 1)
    itr = size(xvals, 2)
    if itr == 0 # x0 is x☆
        x0 = pr.x0
        pDict = pr.p
        f0 = pr.objective(x0, pDict, getGradientToo=false)
        x☆ = x0
        f☆ = f0
    else
        x☆ = xvals[:, itr]
        f☆ = fvals[itr]
    end

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

    fvalMessageMinimum = fvalPrefixMinimum*"$(f☆)"

    myprintln1(v, fvalMessageMinimum, color=:green)

    myprintln1(v, "***************************")
    myprintln1(v, xvalMessage, color=:green)
    if n <= 15
        # If the vector is short enough, print all elements
        for k in 1:n
            myprintln1(v, "x☆[$k] = $(x☆[k])", color=:light_green)
        end
    else
        # If the vector is long, print an even distribution of elements
        myprintln1(v, "Only upto 15 values will be printed.", color=:light_green)

        mid_indices = round.(Int, LinRange(6, n-5, 7)[2:end-1])
        indices_to_print = unique([1:5; mid_indices; n-4:n])
        for k in indices_to_print
            myprintln1(v, "x☆[$k] = $(x☆[k])", color=:light_green)
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

