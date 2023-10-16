using DataFrames

function optimize(pr; 
    verbose::Bool=false, 
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
    fevals += 1
    n = length(x)
    itr = 1
    fvals, αvals = [zeros(Float64, maxiter) for _ in 1:2]
    backtrackVals = zeros(Int64, maxiter, 1)
    xvals = zeros(Float64, n, maxiter)
    
    myprintln(true, "Begin with the solver:", log=log, log_path=log_txt)
    
    while abs(fnext - fₖ) ≥ dftol && itr ≤ maxiter
        printOrNot = verbose && (itr % progress == 0)
        # printOrNot = false
        myprintln(printOrNot, "Iteration $(itr):", log_path=log_txt)
        fₖ, ∇fₖ = obj(x, p)
        fevals += 1
        gevals += 1
        pₖ = findDirection(pr, ∇fₖ)
        
        
        α, x, fnext, backtrackNum, fevals_ls, gevals_ls = (linesearchMethod == "Armijo") ? linesearchArmijo(pr, x, pₖ, itrStart=itrStart, verbose=printOrNot) : strongWolfeBisection(pr, x, pₖ, itrStart=itrStart, verbose=printOrNot)

        myprintln(printOrNot, "Iteration $(itr): x = $(x) is a better point with new fval = $(fnext).", log_path=log_txt)

        fevals += fevals_ls
        gevals += gevals_ls
        fvals[itr] = fnext
        αvals[itr] = α
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
        fvals, αvals, backtrackVals = [arr[1:itr-1] for arr in (fvals, αvals, backtrackVals, xvals)]
        xvals = xvals[:, 1:itr-1]
    end
    
    res = (converged=converged, statusMessage=statusMessage, fvals=fvals, αvals=αvals, backtrackVals=backtrackVals, xvals=xvals, M=M, fevals=fevals, gevals=gevals, pr=pr)

    return res
end


function findDirection(pr::NamedTuple, ∇fnow::Vector{Float64};
    verbose::Bool=false)::Vector{Float64}
    method = pr.alg.method
    n = length(∇fnow)
    if method == "GradientDescent"
        # Bₖ = I(n)
        # pₖ = -Bₖ*∇fnow
        pₖ = -∇fnow
    elseif method == "ConjugateGradientDescent"
        @error "Currently not formulated for this method"
    elseif method == "QuasiNewton"
        @error "Currently not formulated for this method"
    else
        @error "Currently not formulated for this method"
    end

    return pₖ
end

