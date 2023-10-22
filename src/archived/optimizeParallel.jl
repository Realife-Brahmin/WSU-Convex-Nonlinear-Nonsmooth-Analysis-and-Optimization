function optimizeParallel(pr; verbose::Bool=false, log::Bool=true, itrStart::Int64=1)
    # Initial settings
    dftol = pr.alg.dftol
    progress = pr.alg.progress
    maxiter = pr.alg.maxiter
    x0 = pr.x0
    x = x0
    obj = pr.objective
    p = pr.p
    fnext = 1e10
    fₖ = obj(x0, p, getGradientToo=false)
    n = length(x)
    itr = 1
    fvals, αvals = [zeros(Float64, maxiter) for _ in 1:2]
    backtrackVals = zeros(Int64, maxiter, 1)
    xvals = zeros(Float64, n, maxiter)
    
    myprintln(true, "Begin with the solver:")
    
    while abs(fnext - fₖ) ≥ dftol && itr ≤ maxiter
        printOrNot = verbose && (itr % progress == 0)
        # printOrNot = false
        myprintln(printOrNot, "Iteration $(itr):", log=true)
        fₖ, ∇fₖ = obj(x, p)
        pₖ = findDirection(pr, ∇fₖ)
        @warn "Line Search hasn't been parallelized. Only running non-parallelized Line Search."
        α, x, fnext, backtrackNum = linesearch(pr, x, pₖ, itrStart=itrStart)
        # α, x, fnext, backtrackNum = linesearch_parallel(pr, x, pₖ, itrStart=itrStart)
        fvals[itr] = fnext
        αvals[itr] = α
        backtrackVals[itr] = backtrackNum
        xvals[:, itr] = x
        itr += 1
    end
    
    if itr > maxiter
        converged = false
        statusMessage = "Failed to converge despite $(maxiter) iterations! 😢"
        @warn statusMessage
    else
        converged = true
        statusMessage = "Convergence achieved in $(itr) iterations 😄"
        myprintln(true, statusMessage)
        # truncating arrays as they weren't filled to capacity
        fvals, αvals, backtrackVals, xvals = [arr[1:itr] for arr in (fvals, αvals, backtrackVals, xvals)]
    end
    
    res = (converged=converged, statusMessage=statusMessage, fvals=fvals, αvals=αvals, backtrackVals=backtrackVals, xvals=xvals)

    return res
end


