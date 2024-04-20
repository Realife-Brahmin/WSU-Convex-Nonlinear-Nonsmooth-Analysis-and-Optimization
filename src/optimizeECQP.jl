include("projectedGradientConjugateGradient.jl")

function optimizeECQP(pr;
    verbose::Bool=false,
    verbose_ls::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")

    objString = pr.objectiveString
    pECQP = pr.p
    if !haskey(pECQP, :subroutineCall)
        myprintln(verbose, "Calling ECQP Solver Independently.")
        subroutineCall = false
        @unpack G, c, Ae, be = pECQP
    else
        myprintln(verbose, "Calling ECQP Solver as a subroutine for ASQP.")
        @unpack G, c, Ae, be, subroutineCall = pECQP
    end

    @show Ae, be
    verbose = verbose && !subroutineCall

    log_txt = log_path * "log_" * objString * "_" * pr.alg[:method] * "_" * string(pr.alg[:maxiter]) * ".txt"
    if !subroutineCall
        if isfile(log_txt)
            rm(log_txt)
        end # remove logfile if present for the run
    end

    solverState = SolverStateECQPType()

    progress = pr.alg[:progress]
    maxiter = pr.alg[:maxiter]
    etol = pr.alg[:etol]

    x0 = pr.x0
    n = length(x0)
    xk = x0

    fvals = zeros(Float64, maxiter)
    xvals = zeros(Float64, n, maxiter)

    # doing this even though GA requires multiple f evals
    myprintln(verbose, "Starting with initial point x = $(xk).", log=verbose, log_path=log_txt)
    f = pr.objective


    AAT = Ae*transpose(Ae)

    f0 = f(x0, pECQP, getGradientToo=false)
    fk = f0
    solState = SolStatePGCGType(x0, G, c, Ae, fk=f0, etol=etol)

    @unpack fevals = solverState
    fevals += 1
    @pack! solverState = fevals

    myprintln(verbose, "which has fval = $(fk)", log=verbose, log_path=log_txt)

    keepIterationsGoing = true
    causeForStopping = []

    while keepIterationsGoing

        @unpack k = solverState

        printOrNot = verbose && ((k - 1) % progress == 0)
        printOrNot_ECQP = printOrNot & verbose_ls

        myprintln(printOrNot, "Iteration k = $(k)", log=verbose, log_path=log_txt)

        @unpack xk, gk, dk, rk = solState
        # saving the current iterates to solState
        km1, xkm1, gkm1, dkm1, rkm1 = k, xk, gk, dk, rk
        num = transpose(rk)*gk

        if k >= maxiter
            push!(causeForStopping, "Iteration limit reached!")
            keepIterationsGoing = false
            break
        elseif num < etol
            push!(causeForStopping, "Convergence Tolerance Reached")
            keepIterationsGoing = false
            break
        end

        @pack! solState = km1, xkm1, gkm1, dkm1, rkm1

        @unpack actions, fevals = solverState

        xkp1, gkp1, dkp1, rkp1, fevals_1PGCG, actions_1PGCG = solveForNextPGCGIterate(xk, gk, dk, rk, num, G, Ae, AAT, verbose=printOrNot_ECQP) # last two oargs are bogus

        fkp1 = f(xkp1, pECQP, getGradientToo=false)
        fevals += 1
        fevals += fevals_1PGCG # bogus, does nothing
        actions = merge(+, actions, actions_1PGCG) # bogus, does nothing

        @pack! solverState = actions, fevals

        # I prefer to only number a completed iteration, as opposed to numbering an in-process/about-to-begin iteration
        k += 1

        xvals[:, k] = xkp1
        fvals[k] = fkp1 # also incorrect

        xk, gk, dk, rk = xkp1, gkp1, dkp1, rkp1

        @pack! solState = xk, gk, dk, rk

        @pack! solState = k
        @pack! solverState = k

    end

    @unpack k = solverState

    if k ≥ maxiter
        converged = false
        statusMessage = "Failed to converge despite $(maxiter) iterations! 😢"
        myprintln(!subroutineCall, statusMessage, log=log, log_path=log_txt)
        @warn statusMessage
    else
        converged = true
        statusMessage = "Convergence achieved in $(k) iterations 😄"
        myprintln(!subroutineCall, statusMessage, log=log, log_path=log_txt)
    end

    @unpack fevals = solverState

    res = (converged=converged, statusMessage=statusMessage, xvals=xvals, fvals=fvals, fevals=fevals, cause=causeForStopping, pr=pr, solState=solState,
        solverState=solverState)

    res = trim_array(res, k)

    return res

end