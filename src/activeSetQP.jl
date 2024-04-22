include("objfuns/objective.jl")
include("objfuns/functionLP.jl")
include("objfuns/functionQP.jl")
include("optimizeECQP.jl")
include("solveLP.jl")

"""
    getECQPStep(prASQP, solStateASQP; verbose=false, verbose_ls=false)

Calculate a single ECQP step within an ASQP framework to adjust the search direction based on active constraints.

# Arguments
- `prASQP`: The problem record for the ASQP solver, containing the current problem settings and state.
- `solStateASQP`: The current solver state of the ASQP, including variables such as `xk` (current point) and `Wk` (working set of constraints).

# Keyword Arguments
- `verbose::Bool=false`: Enables detailed output throughout the process for debugging and insights.
- `verbose_ls::Bool=false`: Activates detailed logging specifically for line search processes within the ECQP subroutine.

# Description
This function orchestrates a single iteration step of the ECQP subroutine as part of a broader ASQP process. It modifies the active set of constraints to incorporate changes in the working set, constructs an ECQP problem based on these constraints, and then invokes the `optimizeECQP` function to solve it. The solution of the ECQP step informs the adjustment of the search direction in the ASQP solver.

The function aggregates constraints based on the active set indices, formulates an ECQP problem using these constraints, and computes a feasible descent direction based on the solution of the ECQP. It returns the difference between the newly computed solution point and the current point as the search direction for the next ASQP iteration.

# Notes
- The function performs a pre-run check to ensure the feasibility of the input `xk` against the modified constraints `Ae*xk = be`.
- Errors in constraint setup or parameter definitions may lead to infeasible solutions or non-convergence within the ECQP solver.
- `verbose` and `verbose_ls` flags should be used judiciously to avoid excessive output, especially in large-scale optimization problems.

# Example
```julia
# Assuming existing ASQP problem setup and solver state
prASQP = generate_pr(...)  # Generate ASQP problem setup
solStateASQP = SolStateASQPType(...)  # Current state of the ASQP solver

# Compute the ECQP step
pk = getECQPStep(prASQP, solStateASQP, verbose=true)

# Utilize the computed step in the ASQP solver
# Continue with ASQP iterations
```
"""
function getECQPStep(prASQP, solStateASQP;
    verbose::Bool=false,
    verbose_ls::Bool=false,
)

    @unpack etol, itol = pr.alg
    @unpack objective, p, objectiveString = prASQP
    objectiveStringECQP = objectiveString*"_ECQPsubroutine"
    problemTypeECQP = "ECQP"
    methodECQP = "ProjectedGradientCG"

    pDictASQP = p
    @unpack G, c, mE, Ae, be, A, b = pDictASQP
    @unpack xk, Wk = solStateASQP

    myprintln(verbose, "We start for a good feasible descent direction from current point xk = $(xk)")
    WIk = Wk[mE+1:end]
    Ae = vcat(Ae, A[WIk .- mE, :])
    be = vcat(be, b[WIk .- mE])

    # myprintln(verbose, "Let's look at the Ae and be being inserted into ECQP solver in order to represent our current Wk:")
    # @show Ae, be

    subroutineCall = true
    # @show subroutineCall
    pDictECQP = Dict(:G=>G, :c=>c, :Ae=>Ae, :be=>be, :subroutineCall=>subroutineCall)
    # pDictECQP = @packDict "{G, c, Ae, be}"
    # pDictECQP[:subroutineCall] = subroutineCall
    # @show pDictECQP
    # @show pDictECQP
    prECQP = generate_pr(objective, xk, problemType=problemTypeECQP, method=methodECQP, params=pDictECQP, objectiveString=objectiveStringECQP, verbose=false)

    # @show prECQP.p
    res = optimizeECQP(prECQP, verbose=verbose, verbose_ls=verbose_ls, log=false)
    
    xvals = res[:xvals]
    itr = size(xvals, 2)
    if itr == 0 # x0 is xkp1
        xkp1 = xk
    else
        xkp1 = xvals[:, itr]
    end

    pk = xkp1 - xk
    return pk
    
end

"""
    computeLagrangianMultipliersQP(xk, G, c, Awk)

Calculate the Lagrangian multipliers for a given iteration in a quadratic programming problem.

# Arguments
- `xk::Vector{Float64}`: The current point in the parameter space where the Lagrangian multipliers are calculated.
- `G::Matrix{Float64}`: The Hessian matrix of the quadratic part of the objective function.
- `c::Vector{Float64}`: The vector of linear coefficients of the objective function.
- `Awk::Matrix{Float64}`: The active constraint matrix at the current iteration.

# Returns
- `lambdask::Vector{Float64}`: The vector of Lagrangian multipliers corresponding to the active constraints at point `xk`.

# Description
This function computes the Lagrangian multipliers, which provide insights into the sensitivity of the objective function's value at the optimum with respect to constraints. The multipliers are calculated by solving the system of equations formed by the KKT conditions for the active constraints. Specifically, it solves `Awk' * lambda = G*xk + c`, where `lambda` are the Lagrangian multipliers.

# Usage
This function is typically used within an iterative quadratic programming solver where active constraints are being managed. The function helps in assessing which constraints are binding and affecting the objective function.

# Example
```julia
# Example of a small quadratic programming problem
G = [2 0; 0 2]   # Hessian matrix of the objective function
c = [-1; -1]     # Linear coefficients
xk = [0.5; 0.5]  # Current solution estimate
Awk = [1 0; 0 1] # Active constraints (e.g., simple bounds x >= 0, y >= 0)

# Compute the Lagrangian multipliers
lambdask = computeLagrangianMultipliersQP(xk, G, c, Awk)
println("Lagrangian Multipliers: ", lambdask)
```
"""
function computeLagrangianMultipliersQP(xk, G, c, Awk)
    rhs = G*xk + c
    lambdask = transpose(Awk) \ rhs
    return lambdask
end
# c = [1, 3, 5, 2]
# n = length(c)
# x0 = rand(n)
# Ae = [1 1 9 5; 3 5 0 8; 2 0 6 13]
# be = [7, 3, 5]
# mE = length(be)
# A = [0 0 -1 -1]
# b = [-1]
# mI = length(b)
# pLP = Dict(:c => c, :Ae => Ae, :be => be, :mE => mE, :A => A, :b => b, :mI => mI)
# w0 = x0

# objective = LPObjectiveFunction
# # objectiveOriginal = LPObjectiveFunction
# # objectiveString = string(objectiveOriginal)
# objectiveString = "lpTestFunction1"
# params = pLP

# pr = generate_pr(objective, w0, params=params, problemType="LP"; objectiveString=objectiveString)

# xf = computeFeasiblePointForLinearConstraints(pr.p)

# xopt, fopt = solveLP(pr)
