include("LineSearchAlgos.jl")
include("findDirection.jl")

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
    fₖ = obj(x0, p, getGradientToo=false)
    if pr.alg.method == "QuasiNewton"
        QNargs = constructorQNargs(pr, fk=fₖ)
    elseif pr.alg.method == "ConjugateGradientDescent"
        CGargs = constructorCGargs(pr)
    end
    fevals += 1
    n = length(x)
    itr = 1
    fvals, αvals, gmagvals = [zeros(Float64, maxiter) for _ in 1:3]
    backtrackVals = zeros(Int64, maxiter, 1)
    xvals, gvals = [zeros(Float64, n, maxiter) for _ in 1:2]
    
    myprintln(true, "Begin with the solver:", log=log, log_path=log_txt)
    keepIterationsGoing = true
    causeForStopping = []

    while keepIterationsGoing

        printOrNot = verbose && (itr % progress == 0)
        printOrNot_ls = printOrNot & verbose_ls

        myprintln(printOrNot, "Iteration $(itr):", log_path=log_txt)

        fₖ, ∇fₖ = obj(x, p)
        @checkForNaN fₖ
        @checkForNaN ∇fₖ
        gmagval = sum(abs.(∇fₖ))
        fevals += 1
        gevals += 1
        if pr.alg.method == "QuasiNewton"
            QNargs.k = itr
            QNargs.xkp1 = x
            QNargs.fk = fₖ
            QNargs.gkp1 = ∇fₖ
            pₖ, QNargs = findDirection(pr, ∇fₖ, QNargs=QNargs)

        elseif pr.alg.method == "ConjugateGradientDescent"
            CGargs.k = itr
            CGargs.xkp1 = x
            CGargs.gkp1 = ∇fₖ
            pₖ, CGargs = findDirection(pr, ∇fₖ, CGargs=CGargs)

        else
            pₖ = findDirection(pr, ∇fₖ)

        end
        
        α, x, fnext, backtrackNum, fevals_ls, gevals_ls = (linesearchMethod == "Armijo") ? ArmijoBackracking(pr, x, pₖ, itrStart=itrStart, verbose=printOrNot_ls) : StrongWolfeBisection(pr, x, pₖ, itrStart=itrStart, verbose=printOrNot_ls)

        myprintln(printOrNot, "Iteration $(itr): x = $(x) is a better point with new fval = $(fnext).", log_path=log_txt)

        if abs(fnext - fₖ) < dftol
            push!(causeForStopping, "Barely changing fval")
            keepIterationsGoing = false
        end
        if gmagval < gtol
            push!(causeForStopping, "Too small gradient")
            keepIterationsGoing = false
        end
        if itr == maxiter
            push!(causeForStopping, "Too many iterations")
            keepIterationsGoing = false
        end

        fevals += fevals_ls
        gevals += gevals_ls
        fvals[itr] = fnext
        αvals[itr] = α
        gvals[:, itr] = ∇fₖ
        gmagvals[itr] = gmagval
        backtrackVals[itr] = backtrackNum
        xvals[:, itr] = x
        itr += 1
    end
    
    if itr > maxiter
        converged = false
        statusMessage = "Failed to converge despite $(maxiter) iterations! 😢"
        myprintln(true, statusMessage, log=log,  log_path=log_txt)
        @warn statusMessage
    else
        converged = true
        statusMessage = "Convergence achieved in $(itr) iterations 😄"
        myprintln(true, statusMessage, log=log, log_path=log_txt)
        # truncating arrays as they weren't filled to capacity
        fvals, gmagvals, αvals, backtrackVals = [arr[1:itr-1] for arr in (fvals, αvals, backtrackVals, xvals)]
        # xvals = xvals[:, 1:itr-1]
        xvals, gvals = [arr[:, 1:itr-1] for arr in (xvals, gvals)]
    end
    
    res = (converged=converged, statusMessage=statusMessage, fvals=fvals, αvals=αvals, backtrackVals=backtrackVals, xvals=xvals, gmagvals=gmagvals, gvals=gvals, M=M, fevals=fevals, gevals=gevals, cause=causeForStopping, pr=pr)

    return res
end

