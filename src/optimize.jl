include("linesearches.jl")
include("findDirection.jl")
include("types.jl")
include("optimizeASQP.jl")
include("optimizeECQP.jl")
include("optimizeGA.jl")
include("optimizeNM.jl")

function extractBestResults(pr, itr, xvals, fvals)
    # Initialize the optimal solution and function value
    xopt, fopt = nothing, nothing

    # Check if there were any iterations
    if itr == 0
        # If no iterations were done, use the initial values as the optimal values
    x0 = pr.x0
        pDict = pr.p
        f0 = pr.objective(x0, pDict, getGradientToo=false)
        xopt = x0
        fopt = f0
    else
        # If iterations were done, extract the optimal values from the last iteration
        xopt = xvals[:, itr]
        fopt = fvals[itr]
    end

    return xopt, fopt
end

    fk = obj(x0, p, getGradientToo=false)
    myprintln(verbose, "which has fval = $(fk)", log_path=log_txt)
    @pack! solState = fk

    if pr.alg[:method] == "QuasiNewton"
        QNState = QNStateType()
    elseif pr.alg[:method] == "ConjugateGradientDescent"
        CGState = CGStateType()
    elseif pr.alg[:method] == "TrustRegion"
        SR1params = SR1paramsType()
        TRparams = TRparamsType()
    end
    
    fevals += 1
    @pack! solverState = fevals

    n = length(xk)

    fvals, αvals, gmagvals = [zeros(Float64, maxiter) for _ in 1:3]
    backtrackVals = zeros(Int64, maxiter)
    xvals, gvals = [zeros(Float64, n, maxiter) for _ in 1:2]
    
    myprintln(true, "Begin with the solver:", log=log, log_path=log_txt)
    keepIterationsGoing = true
    causeForStopping = []

    justRestarted = false # automatically false if not doing CGD, and if doing CGD and latest β was not zero.

    while keepIterationsGoing

        @unpack k = solverState

        printOrNot = verbose && ( (k - 1) % progress == 0)
        printOrNot_ls = printOrNot & verbose_ls


        myprintln(printOrNot, "Iteration $(k):", log_path=log_txt)

        fk, gk = obj(xk, p)
        @checkForNaN fk
        @checkForNaN gk

        gmagk = sum(abs.(gk))
        usingCGD = false
        fevals += 1
        gevals += 1

        @pack! solState = fk, gk, gmagk
        @pack! solverState = fevals, gevals

        if pr.alg[:method] == "QuasiNewton"
            @pack! QNState = k, xk, fk, gk
            pk, QNState = findDirection(pr, gk, QNState=QNState)

        elseif pr.alg[:method] == "ConjugateGradientDescent"
            usingCGD = true
            @pack! CGState = k, gk
            pk, CGState = findDirection(pr, gk, CGState=CGState)

            @unpack justRestarted = CGState 
        else
            pk = findDirection(pr, gk)
        end
        
        @pack! solState = pk 

        if linesearchMethod == "StrongWolfe"

            solState, solverState = StrongWolfe(pr, solState, solverState,
            verbose=printOrNot_ls)


        elseif linesearchMethod == "Armijo"
            @error "Armijo no longer supported."
        
        else
            @error "Unknown linesearch method"
        end

        @unpack success_ls = solverState
        if ~success_ls
            myprintln(true, "Line search failed... Bad direction or optimal point?")
            push!(causeForStopping, "LineSearch failed.")
            keepIterationsGoing = false
        end

        @unpack xkm1, xk, fkm1, fk, gkm1, gk, gmagkm1, gmagk = solState

        myprintln(printOrNot, "Iteration $(k): x = $(xk) is a better point with new fval = $(fk).", log_path=log_txt)

        # if !usingCGD && !justRestarted && abs(fk - fkm1) < dftol
        #     push!(causeForStopping, "Barely changing fval")
        #     keepIterationsGoing = false
        # end
        if !usingCGD && !justRestarted && gmagkm1 < gtol
            push!(causeForStopping, "Too small gradient at previous step.")
            keepIterationsGoing = false
        end
        if !justRestarted && gmagk < gtol
            push!(causeForStopping, "Too small gradient at latest step.")
            keepIterationsGoing = false
        end
        if k == maxiter
            push!(causeForStopping, "Too many iterations")
            keepIterationsGoing = false
        end

        @unpack Hk, alphak = solState
        @unpack alpha_evals = solverState

        fvals[k] = fk
        αvals[k] = alphak
        gvals[:, k] = gk
        gmagvals[k] = gmagk
        backtrackVals[k] = alpha_evals
        xvals[:, k] = xk

        k += 1

        @pack! solverState = k
        @pack! solState = k

    end
    
    @unpack k = solverState

    if k ≥ maxiter
        converged = false
        statusMessage = "Failed to converge despite $(maxiter) iterations! 😢"
        myprintln(true, statusMessage, log=log,  log_path=log_txt)
        @warn statusMessage
    else
        converged = true
        statusMessage = "Convergence achieved in $(k) iterations 😄"
        myprintln(true, statusMessage, log=log, log_path=log_txt)
    end
    
    res = (converged=converged, statusMessage=statusMessage, fvals=fvals, 
    αvals=αvals, backtrackVals=backtrackVals, xvals=xvals, gmagvals=gmagvals, 
    gvals=gvals, M=M, fevals=fevals, gevals=gevals, cause=causeForStopping, 
    pr=pr)

    res = trim_array(res, k-1)
    return res
end

"""
    trailoredOptimize(pr; nStart=4, factor=2, verbose=false, verbose_ls=false)

Dynamically orchestrates optimization based on the provided problem record and method, adapting parameters and strategies accordingly.

# Arguments
- `pr`: The problem record containing all necessary settings, including the initial point, problem data, algorithm configuration, and optimization method.

# Keyword Arguments
- `nStart::Int=4`: Starting size of the problem dimension for methods requiring warm starts.
- `factor::Int=2`: The factor by which the problem dimension is increased on each warm start iteration.
- `verbose::Bool=false`: Enables verbose output, providing detailed logs of the optimization process.
- `verbose_ls::Bool=false`: Enables verbose output specifically for line search operations.

# Description
This function serves as a universal gateway to various optimization methods, deciding and executing an appropriate strategy based on the method specified in `pr`. It supports a range of methods, from gradient-based to evolutionary algorithms, and includes specific logic for handling warm starts in a scalable manner for certain problem types.

# Notes
- The function is particularly robust in scenarios where adaptive strategies can leverage problem structure and solver capabilities to enhance optimization performance.
- Ensure all required data fields are correctly populated in `pr` to avoid runtime errors.
- Use verbose flags to monitor detailed progression for troubleshooting and performance analysis.
- Each optimization method within the function has specific requirements and behaviors; users should be familiar with these before executing the function.

# Usage
This function is designed to be versatile, catering to different needs and providing flexibility in handling various optimization scenarios, making it especially useful for experimental and complex optimization landscapes.

# Example
```julia
# Define the problem with specific parameters and method
pr = {
    objectiveString: "QuadraticObjective",
    p: {G: Matrix, c: Vector, Ae: Matrix, be: Vector},
    alg: {method: "ActiveSetQP", maxiter: 100, progress: 10},
    x0: [initial_guess]
}

# Execute the adaptive optimization
result = trailoredOptimize(pr, verbose=true)

# Process and utilize the results
println("Optimization result: ", result)
```
"""
function tailoredOptimize(pr; 
    nStart::Int = 4,
    factor::Int = 2,
    verbose=false, 
    verbose_ls=false)

    if warmStart && string(pr.objective) == "drag"
        if pr.alg[:method] == "QuasiNewton"
            nMax = 1024
        elseif pr.alg[:method] == "ConjugateGradientDescent"
            nMax = 1024
        elseif pr.alg[:method] == "GradientDescent"
            nMax = 1024
        else
            @error "Unkown method"
        end

        n = nStart
        x0 = Float64.(collect(LinRange(0.0, 1.0, n+2)[2:n+1]))
        while n <= nMax
            pr = replace_field(pr, :x0, x0)
            # res = optimize(pr, verbose=verbose, verbose_ls=verbose_ls)
            res = optimize2(pr, verbose=verbose, verbose_ls=verbose_ls)

            xopt = res.xvals[:, end]
            x0 = extrapolate(xopt, factor)
            n = length(x0)
        end
    
    elseif pr.alg[:method] == "ActiveSetQP"

        if functionName == "transmissionLines"
            myprintln(verbose, "We're gonna iterate over each area, get the optimal transmission tower locations, and plot them into a figure.", color=:magenta)
            @unpack numOuterAreas = pr.p
            myprintln(verbose, "There are $(numOuterAreas) communities which need to be connected to the central community with the power station.", color=:magenta)
            res = Vector{Any}(undef, numOuterAreas)
            for poly ∈ 1:numOuterAreas
                prPoly = preparePolyProblem(poly=poly)
                @time resPoly = optimizeASQP(prPoly, verbose=verbose, verbose_ls=verbose_ls)
                res[poly] = resPoly
            end
            # @error("Figure out what to do.")
        else
            @time res = optimizeASQP(pr, verbose=verbose, verbose_ls=verbose_ls)
        end

    elseif pr.alg[:method] == "ProjectedGradientCG"
        @time res = optimizeECQP(pr, verbose=verbose, verbose_ls=verbose_ls)
        
    elseif pr.alg[:method] == "NelderMead"
        @time res = optimizeNM(pr, verbose=verbose, verbose_ls=verbose_ls)

    elseif pr.alg[:method] == "GeneticAlgorithm"
        @time res = optimizeGA(pr, verbose=verbose,
        verbose_ls=verbose_ls)
    else
        # @time res = optimize(pr, verbose=verbose, verbose_ls=verbose_ls)
        @time res = optimize2(pr, verbose=verbose, verbose_ls=verbose_ls)
    end
    return res
end

