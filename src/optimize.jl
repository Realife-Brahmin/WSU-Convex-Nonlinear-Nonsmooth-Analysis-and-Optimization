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

"""
    findDirection(pr::NamedTuple, ∇fnow::Vector{Float64}; verbose::Bool=false) -> Vector{Float64}

Compute the search direction for optimization methods based on the provided gradient `∇fnow` and the method specified in `pr.alg.method`.

# Arguments
- `pr::NamedTuple`: A named tuple containing problem configurations. Specifically, it must have `pr.alg.method` which defines the optimization method to be used.
- `∇fnow::Vector{Float64}`: The current gradient of the function to be optimized.

# Keyword Arguments
- `verbose::Bool=false`: Enables additional print statements for debugging and information purposes.

# Returns
- `Vector{Float64}`: The computed search direction.

# Example
```julia
pr = (alg=(method="GradientDescent", ...), ...)
gradient = [1.0, 2.0, 3.0]
direction = findDirection(pr, gradient)
"""
function findDirection(pr::NamedTuple, ∇fnow::Vector{Float64};
    verbose::Bool=false)::Vector{Float64}
    method = pr.alg.method
    n = length(∇fnow)
    if method == "GradientDescent"
        # Bₖ = I(n)
        # pₖ = -Bₖ*∇fnow
        pₖ = -∇fnow
    else 
        @error "Currently not formulated for this method"
    end

    return pₖ
end

"""
    linesearch(pr::NamedTuple, xnow::Vector{Float64}, pₖ::Vector{Float64}; verbose::Bool=false)::Tuple{Float64, Vector{Float64}, Float64}

Performs line search to find an appropriate step size (`α`) that ensures the next parameter value `xnext` satisfies the specified conditions, and returns the objective function value `F` at that point.

# Arguments
- `pr::NamedTuple`: An object containing configurations, data, and algorithm settings.
- `xnow::Vector{Float64}`: The current values of the model parameters.
- `pₖ::Vector{Float64}`: The direction vector for the search.

# Keyword Arguments
- `verbose::Bool`: A flag for printing additional information during execution. Default is `false`.

# Returns
- A tuple containing:
    - `α`: The calculated step size.
    - `x`: The next parameter value `xnow + α*pₖ`.
    - `F`: The objective function value at `x`.

### Notes:
- The specific line search condition to use (e.g., "Armijo" or "StrongWolfe") is specified within the `pr` named tuple.
- This function primarily uses the Armijo condition to determine the step size. 
- It makes use of the `evaluateFunction` and `computeCost` functions.

# Example
```julia
pr = (alg=(linesearch="Armijo", c1=0.1, c2=0.9, ...), ...)
x_values = [1.0, 2.0, 3.0]
direction = [-0.5, -0.5, -0.5]
result = linesearch(pr, x_values, direction, verbose=false)
"""
function linesearch(pr::NamedTuple, xnow::Vector{Float64}, 
    pₖ::Vector{Float64};
    itrMax::Int64=50,
    itrStart::Int64=1,
    verbose::Bool=false,
    log::Bool=true)
    # f = Symbol(pr.objective)
    
    linesearch = pr.alg.linesearch
    c₁ = pr.alg.c1
    c₂ = pr.alg.c2
    β = 1/2^(itrStart-1)
    diff = β*pₖ
    xnext = xnow+diff
    fₖ, ∇fₖ = computeCost(pr, xnow, verbose=verbose, log=log)
    fnext = fₖ
    itr_search_for_α = itrStart-1
    myprintln(verbose, "Current value of F, fₖ = $(fₖ)", log=log)
    armijoSatisfied = false
    strongWolfeSatisfied = false
    if linesearch == "StrongWolfe"
        while !strongWolfeSatisfied && itr_search_for_α ≤ itrMax
            diff = β*pₖ
            myprintln(false, "Let's shift x by $(diff)", log=log)
            xnext = xnow+diff
            fnext = computeCost(pr, xnext, getGradientToo=false)
            # println(c₁*β*∇fₖ'*pₖ)
            myprintln(false, "To be compared against: $(fₖ + c₁*β*∇fₖ'*pₖ)", log=log)
            if fnext ≤ fₖ + c₁*β*∇fₖ'*pₖ
                myprintln(verbose, "Armijo condition satisfied for β = $(β)", log=log)
                fnext, ∇fnext = computeCost(pr, xnext)
                if abs(∇fnext'*pₖ) ≥ abs(c₂*∇fₖ'*pₖ)
                    myprintln(verbose, "Curvature condition satisfied for β = $(β)", log=log)
                    strongWolfeSatisfied = true
                else
                    itr_search_for_α += 1
                    myprintln(false, "Curvature condition NOT satisfied for β = $(β)", log=log)
                    β /= 2
                    myprintln(false, "Line Search Iterations = $(itr_search_for_α)", log=log)
                end
            else
                itr_search_for_α += 1
                myprintln(verbose, "Armijo condition NOT satisfied for β = $(β)", log=log)
                β /= 2
                myprintln(verbose, "Line Search Iterations = $(itr_search_for_α)", log=log)
            end 
        end
    elseif linesearch == "Armijo"
        # fₖ, ∇fₖ = pr.objective( xnow, t)
        while !armijoSatisfied && itr_search_for_α ≤ itrMax
            diff = β*pₖ
            myprintln(verbose, "Let's shift x by $(diff)", log=log)
            xnext = xnow+diff
            fnext = computeCost(pr, xnext, getGradientToo=false)
            # println(c₁*β*∇fₖ'*pₖ)
            myprintln(verbose, "To be compared against: $(fₖ + c₁*β*∇fₖ'*pₖ)", log=log)
            if fnext ≤ fₖ + c₁*β*∇fₖ'*pₖ
                myprintln(verbose, "Armijo condition satisfied for β = $(β)", log=log)
                armijoSatisfied = true
            else
                itr_search_for_α += 1
                myprintln(verbose, "Armijo condition NOT satisfied for β = $(β)", log=log)
                β /= 2
                myprintln(verbose, "Line Search Iterations = $(itr_search_for_α)", log=log)
            end 
        end
    else 
        @error "Unknown linesearch condition"
    end
    
    if itr_search_for_α > itrMax
        @error "Line Search failed at point x = $(xnext) despite $(itr_search_for_α) iterations."
    end

    α = β
    return (α=α, x=xnext, f=fnext, backtracks=itr_search_for_α) 
end

# end