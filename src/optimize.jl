include("LineSearchAlgos.jl")
include("linesearches.jl")
include("findDirection.jl")
include("types.jl")

function optimize(pr; 
    verbose::Bool=false, 
    verbose_ls::Bool=false,
    log::Bool=true,
    log_path::String="./logging/",
    itrStart::Int64=1)

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
    x = x0

    myprintln(verbose, "Starting with initial point x = $(x).", log_path=log_txt)
    obj = pr.objective
    p = pr.p
    M = max(size(p.data, 1), 1)
    fnext = 1e10
    fk = obj(x0, p, getGradientToo=false)
    solState.fk = fk
    if pr.alg.method == "QuasiNewton"
        QNargs = constructorQNargs(pr, fk=fk)
    elseif pr.alg.method == "ConjugateGradientDescent"
        CGargs = constructorCGargs(pr)
    end
    
    fevals += 1
    @pack! solverState = fevals
    n = length(x)

    fvals, αvals, gmagvals = [zeros(Float64, maxiter) for _ in 1:3]
    backtrackVals = zeros(Int64, maxiter)
    xvals, gvals = [zeros(Float64, n, maxiter) for _ in 1:2]
    
    myprintln(true, "Begin with the solver:", log=log, log_path=log_txt)
    keepIterationsGoing = true
    causeForStopping = []

    CGDRestartFlag = false # automatically false if not doing CGD, and if doing CGD and latest β was not zero.

    while keepIterationsGoing

        k = solverState.k

        printOrNot = verbose && (k % progress == 0)
        printOrNot_ls = printOrNot & verbose_ls


        myprintln(printOrNot, "Iteration $(k):", log_path=log_txt)

        fk, gk = obj(x, p)
        @checkForNaN fk
        @checkForNaN gk

        gmagk = sum(abs.(gk))
        
        fevals += 1
        gevals += 1

        @pack! solState = fk, gk, gmagk
        @pack! solverState = fevals, gevals

        if pr.alg.method == "QuasiNewton"
            QNargs.k = k
            QNargs.xkp1 = x
            QNargs.fk = fk
            QNargs.gkp1 = gk
            pk, QNargs = findDirection(pr, gk, QNargs=QNargs)

        elseif pr.alg.method == "ConjugateGradientDescent"
            CGargs.k = k
            pk, CGargs = findDirection(pr, gk, CGargs=CGargs)
            # CGDRestartFlag = CGargs.justRestarted
            CGDRestartFlag = false # temporary until new types are inserted
        else
            pk = findDirection(pr, gk)

        end
        
        @pack! solState = pk 

        if linesearchMethod == "StrongWolfe"

            solState, solverState = StrongWolfe(pr, solState, solverState, verbose=printOrNot_ls) # under construction
            # if success == false && pr.alg.method == "ConjugateGradientDescent"
            #     CGDRestartFlag = true
            # end

        elseif linesearchMethod == "Armijo"
            @error "Armijo no longer supported."
        
        else
            @error "Unknown linesearch method"
        end

        @unpack xkm1, xk, fkm1, fk, gkm1, gk, gmagkm1, gmagk = solState

        myprintln(printOrNot, "Iteration $(k): x = $(x) is a better point with new fval = $(fk).", log_path=log_txt)

        if !CGDRestartFlag && abs(fk - fkm1) < dftol
            push!(causeForStopping, "Barely changing fval")
            keepIterationsGoing = false
        end
        if !CGDRestartFlag && gmagkm1 < gtol
            push!(causeForStopping, "Too small gradient at previous step.")
            keepIterationsGoing = false
        end
        if !CGDRestartFlag && gmagk < gtol
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
    
    res = (converged=converged, statusMessage=statusMessage, fvals=fvals, αvals=αvals, backtrackVals=backtrackVals, xvals=xvals, gmagvals=gmagvals, gvals=gvals, M=M, fevals=fevals, gevals=gevals, cause=causeForStopping, pr=pr)

    res = trim_array(res, k-1)
    return res
end

