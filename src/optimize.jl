# module optimize

# export optimize

function optimize(pr; verbose::Bool=false, log::Bool=true, itrStart::Int64=1)
    # Initial settings
    dftol = pr.alg.dftol
    progress = pr.alg.progress
    maxiter = pr.alg.maxiter
    x0 = pr.x0
    x = x0
    fnext = 1e10
    fₖ = computeCost(pr, x0, getGradientToo=false)
    n = length(x)
    itr = 1
    fvals, αvals = [zeros(Float64, maxiter) for _ in 1:2]
    backtrackVals = zeros(Int64, maxiter, 1)
    xvals = zeros(Float64, n, maxiter)
    
    myprintln(verbose, "Begin with the solver:")
    
    while abs(fnext - fₖ) ≥ dftol && itr ≤ maxiter
        printOrNot = verbose && (itr % progress == 0)
        myprintln(printOrNot, "Iteration $(itr):", log=true)
        fₖ, ∇fₖ = computeCost(pr, x)
        pₖ = findDirection(pr, ∇fₖ)
        α, x, fnext, backtrackNum = linesearch(pr, x, pₖ, verbose=printOrNot, itrStart=itrStart)
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
        println(statusMessage)
        # truncating arrays as they weren't filled to capacity
        fvals, αvals, backtrackVals, xvals = [arr[1:itr] for arr in (fvals, αvals, backtrackVals, xvals)]
    end
    
    res = (converged=converged, statusMessage=statusMessage, fvals=fvals, αvals=αvals, backtrackVals=backtrackVals, xvals=xvals)

    return res
end

# end