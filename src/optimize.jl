include("linesearches.jl")
include("findDirection.jl")
include("types.jl")

function optimize(pr; 
    verbose::Bool=false, 
    verbose_ls::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")

    log_txt = log_path*"log_"*string(pr.objective)*"_"*pr.alg.method*"_"*pr.alg.linesearch*"_"*string(pr.alg.maxiter)*".txt"

    if isfile(log_txt)
        rm(log_txt)
    end # remove logfile if present for the run

    solverState = SolverStateType()
    solState = SolStateType(xk=pr.x0)

    # Initial settings
    fevals = 0
    gevals = 0
    dftol = pr.alg.dftol
    gtol = pr.alg.gtol
    progress = pr.alg.progress
    maxiter = pr.alg.maxiter
    linesearchMethod = pr.alg.linesearch
    x0 = pr.x0
    xk = x0

    myprintln(verbose, "Starting with initial point x = $(xk).", log_path=log_txt)
    obj = pr.objective
    p = pr.p
    data = p[:data]
    M = max(size(data, 1), 1)

    fk = obj(x0, p, getGradientToo=false)
    myprintln(verbose, "which has fval = $(fk)", log_path=log_txt)
    @pack! solState = fk

    if pr.alg.method == "QuasiNewton"
        QNState = QNStateType()
    elseif pr.alg.method == "ConjugateGradientDescent"
        CGState = CGStateType()
    elseif pr.alg.method == "TrustRegion"
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

        if pr.alg.method == "QuasiNewton"
            @pack! QNState = k, xk, fk, gk
            pk, QNState = findDirection(pr, gk, QNState=QNState)

        elseif pr.alg.method == "ConjugateGradientDescent"
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

        if !usingCGD && !justRestarted && abs(fk - fkm1) < dftol
            push!(causeForStopping, "Barely changing fval")
            keepIterationsGoing = false
        end
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
    warm_start_optimize(pr; nStart::Int = 4, factor::Int = 2, verbose=false, verbose_ls=false) -> NamedTuple

Conducts warm start optimization for estimating a monotonic function that minimizes a drag function, progressively refining the solution by increasing the number of points defining the function. This approach is currently specialized for the 'drag' objective function and will perform a standard single-shot optimization for other objectives.

# Arguments
- `pr`: A NamedTuple containing the problem definition and settings for optimization.
- `nStart::Int`: The starting number of points for the warm start optimization process.
- `factor::Int`: The factor by which the number of points is increased in each extrapolation step.
- `verbose`: A Boolean flag for verbose output during optimization.
- `verbose_ls`: A Boolean flag for verbose output during line search.

# Behavior
For the 'drag' objective function, the method starts with `nStart` points, optimizes, and then uses the result to extrapolate a finer initial guess for the next round, increasing the number of points each time by the specified `factor`. This continues until `nMax` points are reached, which depends on the optimization method used:
- `QuasiNewton`: `nMax` is set to 1024.
- `ConjugateGradientDescent`: `nMax` is set to 1024.
- Other methods default to `nMax` of 512.

For objective functions other than 'drag', `warm_start_optimize` performs a single-shot optimization using the `optimize` function.

# Returns
- `res`: A NamedTuple containing the results of the optimization process.

# Notes
This function is designed for use with optimization problems where a good initial guess can significantly speed up convergence. The current implementation is tuned for the 'drag' function; using it with other objectives will not leverage the warm start capability.

Example usage:
```julia
pr = (objective=drag, alg=alg, x0=[0.25, 0.5, 0.75], ...)
result = warm_start_optimize(pr, verbose=true, verbose_ls=false)
"""
function warm_start_optimize(pr; 
    nStart::Int = 4,
    factor::Int = 2,
    verbose=false, 
    verbose_ls=false)

    if warmStart && string(pr.objective) == "drag"
        if pr.alg.method == "QuasiNewton"
            nMax = 1024
        elseif pr.alg.method == "ConjugateGradientDescent"
            nMax = 1024
        elseif pr.alg.method == "GradientDescent"
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
    else
        @time res = optimize(pr, verbose=verbose, verbose_ls=verbose_ls)
    end
    return res
end

