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
        # α, x, fnext, backtrackNum = linesearch(pr, x, pₖ, itrStart=itrStart)
        α, x, fnext, backtrackNum = linesearch_parallel(pr, x, pₖ, itrStart=itrStart)
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

function linesearch_parallel(pr::NamedTuple, xnow::Vector{Float64}, 
    pₖ::Vector{Float64};
    itrMax::Int64=50,
    itrStart::Int64=1,
    verbose::Bool=false,
    log::Bool=true)

    obj = pr.objective
    p = pr.p
    isStrongWolfe = (pr.alg.linesearch == "StrongWolfe")
    c₁ = pr.alg.c1
    fₖ, ∇fₖ = obj(xnow, p)
    fnext = fₖ

    # Atomic variable to store the result.
    result = Threads.Atomic{Any}(nothing)

    @threads for itr_search_for_α in itrStart:itrMax
        # Check if result is already found by other thread
        if isnothing(result[])
            β = 1 / 2^(itr_search_for_α-1)
            xnext = copy(xnow)
            @inbounds xnext .= xnow .+ β .* pₖ
            fnext, ∇fnext = obj(xnext, p)
            comparison_val = fₖ + c₁ * β * dot(∇fₖ, pₖ)

            if fnext ≤ comparison_val
                if isStrongWolfe && abs(dot(∇fnext, pₖ)) < abs(c₁ * dot(∇fₖ, pₖ))
                    continue
                else
                    # Update the result atomically
                    Threads.atomic_cas!(result, nothing, (β, xnext, fnext, itr_search_for_α))
                end
            end
        end
    end

    # Extract the result from atomic variable
    res = result[]
    if isnothing(res)
        @error "Line Search failed at point x = $(xnow) despite $(itrMax) iterations."
    else
        α, x, f, backtracks = res
        return (α=α, x=x, f=f, backtracks=backtracks)
    end
end