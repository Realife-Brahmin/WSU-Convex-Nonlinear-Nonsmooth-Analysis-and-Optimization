using Parameters

function SolStateALPType(xk; # xk is actually wk = [xk (original variables); yk (inequality slack variables)]
    fkm1=100.0,
    fk=100.0,
    etol=1e-8,
    gtol=1e-8
    )
    # N = length(xk) # actually n+mI
    
    # Prepare the state dictionary with initial values
    solState = Dict(
        :km1 => -1, :k => 0,
        :xkm1 => myfill(xk, -27.0), :xk => xk,
        :fkm1 => fkm1, :fk => fk,
        :gkm1 => myfill(xk, 11.3), :gk => myfill(xk, 22.7),
        :etol => etol,
        :gtol => gtol
    )

    return solState
end

function SolverStateALPType(;
    k=0,
    fevals=0,
    gevals=0,
    lagrevals=0,
    actions=Dict()
)

    solverState = Dict(:k => k, :fevals => fevals, :gevals=>gevals, :lagrevals => lagrevals, :actions => actions)

    return solverState
end

"""
    SolStateASQPType(xk, Ae; Wk0=[], fkm1=100.0, fk=100.0, itol=1e-8)

Initialize the solution state for an Active Set Quadratic Programming (ASQP) solver.

# Arguments
- `xk::Vector{Float64}`: The initial point in the parameter space where the solver starts.
- `Ae::Matrix{Float64}`: The matrix representing equality constraints.

# Keyword Arguments
- `Wk0::Vector{Int}=[]`: Initial set of active (working) constraints indices. If empty, defaults to indices of all equality constraints.
- `fkm1::Float64=100.0`: The function value at the previous iteration (initially set to a high value to represent uninitialized state).
- `fk::Float64=100.0`: The function value at the current iteration (initially set to a high value to represent uninitialized state).
- `itol::Float64=1e-8`: The iteration tolerance for convergence checking.

# Returns
- `solState::Dict`: A dictionary containing the state of the ASQP solver, encapsulating various metrics and values such as iteration counts, function values, and active constraint sets.

# Description
This function sets up and returns a dictionary encapsulating the necessary state for the ASQP algorithm. It initializes the state with default values for function evaluations and active constraints based on the provided equality constraints matrix `Ae`. This state aids in guiding the iterative process of the ASQP solver by tracking previous and current values necessary for convergence checks and updates.

# Usage
The `SolStateASQPType` is particularly useful at the beginning of an optimization process to set up the required state information for managing constraint activation and deactivation across iterations.

# Example
```julia
# Define initial conditions and problem parameters
xk = [0.5, 0.5]  # Initial guess for variables
Ae = [1 0; 0 1]  # Equality constraint coefficients matrix

# Initialize the solver state for an ASQP optimization
solState = SolStateASQPType(xk, Ae)

# Example output of the initial solver state
println("Initial Solver State: ", solState)
```
"""
function SolStateASQPType(xk, Ae;
    Wk0=[],
    fkm1=100.0,
    fk=100.0,
    itol=1e-8)

    mE = size(Ae, 1)
    
    if isempty(Wk0)
        Wk0 = collect(1:mE) # at least all the equality constraints (may be zero => still empty)
    end

    # Prepare the state dictionary with initial values
    solState = Dict(
        :km1 => -1, :k => 0,
        :xkm1 => myfill(xk, -27.0), :xk => xk,
        :fkm1 => fkm1, :fk => fk,
        :Wkm1 => -1 * Wk0, # has no sense in reality
        :Wk => Wk0,
        :itol => itol
    )

    return solState
end

"""
    SolverStateASQPType(; k=0, fevals=0, lagrevals=0, actions=Dict())

Initialize the solver state for an Active Set Quadratic Programming (ASQP) method.

# Keyword Arguments
- `k::Int=0`: The current iteration number.
- `fevals::Int=0`: The number of function evaluations performed.
- `lagrevals::Int=0`: The number of Lagrangian multiplier evaluations performed.
- `actions::Dict=Dict()`: A dictionary to store additional information about actions taken during the solver's execution.

# Returns
- `solverState::Dict`: A dictionary containing the state of the ASQP solver, including iteration counts and evaluations.

# Description
This function creates and returns a dictionary that records various aspects of the state of an ASQP solver, such as the iteration count, number of function evaluations, number of Lagrangian multiplier evaluations, and other actions. This state is used to track progress and manage the execution flow of the optimization process.

# Usage
`SolverStateASQPType` is typically used at the initialization stage of an ASQP solver to set up the initial state of the solver. It can be updated in subsequent iterations to reflect ongoing progress and actions.

# Example
```julia
# Initialize the solver state for an ASQP optimization
solverState = SolverStateASQPType()

# Print the initial state of the solver
println("Initial Solver State: ", solverState)
```
"""
function SolverStateASQPType(;
    k=0,
    fevals=0,
    lagrevals=0,
    actions=Dict()
)

    solverState = Dict(:k => k, :fevals => fevals, :lagrevals => lagrevals, :actions => actions)

    return solverState
end

"""
    SolStateType

A mutable struct used to maintain the state of the solution in an optimization solver.

# Fields
- `k::Int`: The current iteration number.
- `xkm1::Vector{Float64}`: The solution vector from the previous iteration.
- `xk::Vector{Float64}`: The current solution vector.
- `fkm1::Float64`: The objective function value at the previous iteration.
- `fk::Float64`: The current objective function value.
- `gkm1::Vector{Float64}`: The gradient of the objective function at the previous iteration.
- `gk::Vector{Float64}`: The current gradient of the objective function.
- `gmagkm1::Float64`: The magnitude of the gradient at the previous iteration.
- `gmagk::Float64`: The current magnitude of the gradient.
- `pkm1::Vector{Float64}`: The search direction used in the previous iteration.
- `pk::Vector{Float64}`: The current search direction.
- `Hkm1::Matrix{Float64}`: The inverse approximate Hessian from the previous iteration.
- `Hk::Matrix{Float64}`: The current approximate inverse Hessian.
- `alphak::Float64`: The step size used in the current iteration.

# Constructor
The constructor for `SolStateType` can be called with keyword arguments corresponding to each field. If a field is not specified, it defaults to a pre-defined value such as `1` for `k`, `0.0` for scalar values, an empty `Float64` array for vector fields, or an undefined `0x0` `Float64` matrix for matrix fields.

# Example Usage
```julia
    solState = SolStateType(k=5, xkm1=[1.0, 2.0], xk=[1.5, 2.5])
```
This will create a `SolStateType` instance with the iteration number set to `5`, previous solution vector as `[1.0, 2.0]`, and current solution vector as `[1.5, 2.5]`. All other fields will be set to their default values.

Please note that the fields `Hkm1` and `Hk` are initialized with undefined matrices. It is expected that the user will provide properly sized and initialized matrices if these fields are to be used.

"""
mutable struct SolStateType
    k::Int
    xkm1::Vector{Float64}
    xk::Vector{Float64}
    fkm1::Float64
    fk::Float64
    gkm1::Vector{Float64}
    gk::Vector{Float64}
    gmagkm1::Float64
    gmagk::Float64
    pkm1::Vector{Float64}
    pk::Vector{Float64}
    Hkm1::Matrix{Float64}
    Hk::Matrix{Float64}
    alphak::Float64
    

    function SolStateType(;
        k=1, 
        xkm1=Float64[], xk=Float64[],
        fkm1=0.0, fk=0.0,
        gkm1=Float64[], gk=Float64[],
        gmagkm1=0.0, gmagk=0.0, 
        pkm1=Float64[], pk=Float64[], 
        Hkm1=Matrix{Float64}(undef, 0, 0), Hk=Matrix{Float64}(undef, 0, 0), 
        alphak=0.0)

        new(k, xkm1, xk, fkm1, fk, gkm1, gk, gmagkm1, gmagk, pkm1, pk, Hkm1, Hk, alphak)
    end

end


"""
    SolStateNMType(; k=0, Xkm1=zeros(0, 0), Xk=zeros(0, 0), fkm1=100.0, fk=100.0, Deltak=100.0)

Initializes and returns a dictionary representing the state of the solution in the Nelder-Mead optimization process.

# Arguments
- `k`: Iteration counter, indicating the current step of the optimization process.
- `Xkm1`: Simplex matrix from the previous iteration (`k-1`), an `n x (n+1)` matrix where `n` is the dimensionality of the space.
- `Xk`: Current simplex matrix, an `n x (n+1)` matrix.
- `fkm1`: The objective function value associated with the best point in the simplex `Xkm1`.
- `fk`: The objective function value associated with the best point in the current simplex `Xk`.
- `Deltak`: The size of the simplex Xk.

# Returns
- A dictionary containing the current state of the optimization process, with keys corresponding to each argument.

# Examples
```julia
# Initialize a solution state for the Nelder-Mead optimization
solState = SolStateNMType(k=1, Xkm1=rand(2, 3), Xk=rand(2, 3), fkm1=5.0, fk=4.5, Delta=0.5)
```
"""
function SolStateNMType(; k=0, Xkm1=zeros(0, 0), Xk=zeros(0, 0),
    fkm1=100.0, fk=100.0, Deltak=100.0)

    return Dict(:k => k, :Xkm1 => Xkm1, :Xk => Xk, :fkm1 => fkm1, :fk => fk, :Deltak => Deltak)

end

"""
    SolverStateNMType(; k=0, fevals=0, actions=Dict(...)) -> Dict

Create and return a dictionary representing the state of a solver at a given iteration within a numerical method algorithm, specifically tailored for algorithms like the Nelder-Mead method. This function captures both the iteration number, the number of function evaluations, and a detailed account of various actions taken during the algorithm's execution.

# Arguments
- `k::Int=0`: The current iteration number. Defaults to 0.
- `fevals::Int=0`: The total number of function evaluations performed up to this point. Defaults to 0.
- `actions::Dict`: A dictionary mapping action names to the number of times each action was performed. The default actions and their initial counts are:
    - `:extend => 0`: Counts how many times the algorithm performed an extension step.
    - `:insideContract => 0`: Counts the number of inside contraction steps.
    - `:outsideContract => 0`: Counts the number of outside contraction steps.
    - `:reflect => 0`: Counts the number of reflection steps performed.
    - `:shrink => 0`: Counts the number of shrink steps performed.
    - `:sort => 0`: Counts the number of sorting operations performed on the vertices of the simplex. 
    - `:insertIntoSorted => 0`: Counts the number of times a point was added to an already sorted simplex. 

# Returns
- `Dict`: A dictionary with keys `:k`, `:fevals`, and `:actions`, corresponding to the function's arguments.

# Example
```julia
solver_state = SolverStateNMType(k=1, fevals=5, actions=Dict(:extend => 1, :reflect => 2))
```
"""
function SolverStateNMType(; 
    k=0, 
    fevals=0,
    actions=Dict(:extend => 0, :extensionFailure => 0 , :extensionSuccess => 0,
    :insideContract => 0, :insideContractFailure => 0, :insideContractSuccess => 0, :outsideContract => 0, 
    :outsideContractFailure =>0, :outsideContractSuccess =>0,
    :reflect => 0, :shrink => 0, :sort => 0, :insertIntoSorted => 0)
    )

    return Dict(:k => k, :fevals => fevals, :actions => actions)

end

"""
    SolStatePGCGType(xk, G, c, Ae; fkm1=100.0, fk=100.0, etol=1e-8)

Initialize the solution state for the Projected Gradient Conjugate Gradient (PGCG) method.

# Arguments
- `xk::Vector{Float64}`: The initial point or current point in the parameter space.
- `G::Matrix{Float64}`: The Hessian matrix or an approximation of the Hessian of the quadratic part of the objective function.
- `c::Vector{Float64}`: The vector of linear coefficients of the objective function.
- `Ae::Matrix{Float64}`: The matrix representing equality constraints.

# Keyword Arguments
- `fkm1::Float64=100.0`: The function value at the previous iteration.
- `fk::Float64=100.0`: The function value at the current iteration.
- `etol::Float64=1e-8`: The error tolerance for convergence checking.

# Returns
- `solState::Dict`: A dictionary containing the state of the optimization process, including various vectors and parameters like gradients, residuals, and function values.

# Description
This function sets up and returns a dictionary encapsulating the state needed for the PGCG algorithm. It initializes residuals, gradients, and search directions based on the provided matrix and vector inputs. This state aids in guiding subsequent iterations of the PGCG solver by tracking previous and current values necessary for convergence checks and updates.

# Usage
The `SolStatePGCGType` is typically called at the beginning of an optimization process to initialize necessary variables and state information.

# Example
```julia
# Define initial conditions and problem parameters
xk = [0.5, 0.5]  # Initial guess
G = [2 0; 0 2]   # Hessian matrix of the objective function
c = [-1; -1]     # Linear coefficients of the objective function
Ae = [1 0; 0 1]  # Equality constraint coefficients

# Initialize the solver state
solState = SolStatePGCGType(xk, G, c, Ae)

# Example output of the initial solver state
println("Initial Solver State: ", solState)
```
"""
function SolStatePGCGType(xk, G, c, Ae; 
    fkm1=100.0,
    fk=100.0,
    etol=1e-8)

    # note that here the variable 'xk' here might be 'wk', which can have a different meaning than the original problem.

    rk = G * xk + c

    AAT_inv = inv(Ae * transpose(Ae))
    gk = rk - transpose(Ae) * AAT_inv * Ae * rk

    # Initialize dk
    dk = -gk

    # Prepare the state dictionary with initial values
    solState = Dict(
        :km1 => -1, :k => 0,
        :xkm1 => myzeros(xk), :xk => xk,
        :rkm1 => myfill(rk, 22), :rk => rk,
        :gkm1 => myfill(gk, 77), :gk => gk,
        :dkm1 => myfill(dk, 33), :dk => dk,
        :fkm1 => fkm1, :fk => fk,
        :etol => etol
    )

    return solState
end

"""
    SolverStateECQPType(; k=0, fevals=0, actions=Dict()) -> Dict

Initializes the state dictionary for an Extended Constrained Quadratic Programming (ECQP) solver. The state dictionary contains key metrics and action records that define the solver's current status and history of operations.

# Keyword Arguments
- `k::Int=0`: The current iteration index in the ECQP solver. Initializes to `0`, representing the starting iteration.
- `fevals::Int=0`: The total number of objective function evaluations performed up to the current iteration. Useful for monitoring the computational cost of the solving process.
- `actions::Dict`: A dictionary with keys representing different actions or events in the ECQP solving process and values indicating the count of each action/event. This could include actions such as "gradient evaluations", "constraint evaluations", "line searches", etc., depending on what is relevant to the ECQP solver's implementation.

# Returns
- `Dict`: A dictionary object encapsulating the ECQP solver's state, including iteration count, function evaluations, and actions taken.

# Example
```julia
# Initialize an ECQP solver state with default values
solverState = SolverStateECQPType()

# Initialize an ECQP solver state with a custom iteration index and actions
customActions = Dict("gradient evaluations" => 10, "line searches" => 5)
solverStateCustom = SolverStateECQPType(k=3, fevals=15, actions=customActions)
```
"""
function SolverStateECQPType(;
    k=0,
    fevals=0,
    actions=Dict()
    )

    return Dict(:k => k, :fevals => fevals, :actions => actions)

end

"""
    SolStateGAType(; k=0, Xkm1=zeros(0, 0), Xk=zeros(0, 0),
                    Fkm1=zeros(0), Fk=zeros(0), fkm1=100.0, fk=100.0)

Construct a dictionary that encapsulates the solution state for a Genetic Algorithm (GA) at a given generation. This state includes information about the population and fitness scores.

# Arguments
- `k::Int=0`: The current generation number. Defaults to 0 for the initial generation.
- `Xkm1::Matrix{Float64}`: The population matrix from the previous generation (k-1). Each row represents an individual in the population. Defaults to an empty matrix.
- `Xk::Matrix{Float64}`: The current population matrix (at generation k). Analogous to `Xkm1`, each row represents an individual. Defaults to an empty matrix.
- `Fkm1::Vector{Float64}`: The fitness vector for the previous generation's population. Each element corresponds to the fitness of the respective individual in `Xkm1`. Defaults to an empty vector.
- `Fk::Vector{Float64}`: The fitness vector for the current generation's population. Each element corresponds to the fitness of the respective individual in `Xk`. Defaults to an empty vector.
- `fkm1::Float64=100.0`: The best fitness value obtained in the previous generation. Defaults to 100.0.
- `fk::Float64=100.0`: The best fitness value obtained in the current generation. Defaults to 100.0.

# Returns
- `Dict`: A dictionary object with keys corresponding to the argument names and their associated values, representing the state of the GA solver.

# Example
```julia
# Initialize state with default values
solverState = SolStateGAType()

# Initialize state with a specific generation k and corresponding data
solverState = SolStateGAType(k=3, Xkm1=rand(10, 5), Xk=rand(10, 5),
                            Fkm1=rand(10), Fk=rand(10), fkm1=0.5, fk=0.3)
```
"""
function SolStateGAType(; 
    k=0, Xkm1=zeros(0, 0), Xk=zeros(0, 0),
    Fkm1=zeros(0), Fk=zeros(0), fkm1=100.0, fk=100.0
    )

    return Dict(:k => k, :Xkm1 => Xkm1, :Xk => Xk, :Fkm1 => Fkm1, :Fk => Fk, :fkm1 => fkm1, :fk => fk)

end

"""
    SolverStateGAType(; k=0, fevals=0, fvalRepeats=0, actions=Dict())

Create a dictionary representing the state of a Genetic Algorithm (GA) solver. This state captures key metrics and event counts at a specific generation of the GA.

# Arguments
- `k::Int=0`: The current generation index. Initializes to 0, representing the starting generation.
- `fevals::Int=0`: The cumulative number of fitness function evaluations conducted up to this generation.
- `fvalRepeats::Int=0`: The count of fitness evaluations that resulted in the same fitness value, indicating potential stagnation in the population.
- `actions::Dict`: A dictionary with keys representing different actions or events in the GA process and values indicating the count of each action/event. The default dictionary includes the following keys:
    - `:bothParentsSurvived`: Count of occurrences where both parents survived to the next generation.
    - `:crossover`: Count of crossover operations performed.
    - `:genFitnessImproved`: Count of generations where fitness has improved.
    - `:genFitnessNotImproved`: Count of generations where fitness has not improved.
    - `:mutation`: Count of mutation operations performed.
    - `:noParentSurvived`: Count of occurrences where neither parent survived to the next generation.
    - `:parentsSelected`: Count of times parents have been selected for mating.
    - `:onlyOneParentSurvived`: Count of occurrences where only one parent survived to the next generation.

# Returns
- `Dict`: A dictionary object encapsulating the GA solver's state.

# Example
```julia
solverState = SolverStateGAType(k=5, fevals=100, fvalRepeats=10)
```
"""
function SolverStateGAType(;
    k=0,
    fevals=0,
    fvalRepeats=0,
    actions=Dict(
        :bothParentsSurvived => 0, 
        :crossover => 0, 
        :genFitnessImproved => 0,
        :genFitnessNotImproved => 0,
        :mutation => 0,
        :noParentSurvived => 0,
        :parentsSelected => 0,
        :onlyOneParentSurvived => 0
        )
    )

    return Dict(:k => k, :fevals => fevals, :fvalRepeats => fvalRepeats, :actions => actions)

end

"""
    SolverStateType

A mutable struct that holds the state of the optimization solver. It tracks the iteration count and the number of function, gradient, and Hessian evaluations, as well as the success of the line search.

Fields:

- `k`: The current iteration number.
- `fevals`: The cumulative number of function evaluations performed.
- `gevals`: The cumulative number of gradient evaluations performed.
- `Hevals`: The cumulative number of Hessian evaluations performed.
- `alpha_evals`: The number of evaluations for the current step size.
- `success_ls`: A Boolean flag indicating whether the last line search was successful.

Constructor:

The constructor accepts keyword arguments for each of the fields. Default values are provided for all fields.
Example Usage:

Create a new instance by specifying any of the fields. For example:
```julia
solverState = SolverStateType(k=10, fevals=50, gevals=50)
```
"""
mutable struct SolverStateType
    k::Int
    fevals::Int
    gevals::Int
    Hevals::Int
    alpha_evals::Int # only current evals
    success_ls::Bool

    function SolverStateType(;
        k=1, 
        fevals=0, 
        gevals=0, 
        Hevals=0, 
        success_ls=false)
        new(k, fevals, gevals, Hevals, success_ls)
    end
end

"""
    InterpolParams

A mutable struct used within the line search algorithm to manage interpolation parameters, specifically focused on determining an acceptable step size, alphaj.

Fields:

- `j`: The current interpolation iteration.
- `alphaj`: The trial step size being evaluated.
- `alphaLo`: The current lower bound on the step size.
- `alphaHi`: The current upper bound on the step size.
- `alphatol`: The tolerance within which the step size is considered acceptable.
- `alphatolBreached`: A Boolean flag indicating whether the tolerance criterion has been breached.
- `change`: A string indicating if the next step size should be increased, decreased, or remain unchanged.

Constructor:

The constructor for InterpolParams accepts keyword arguments for each of the fields. Default values are assumed for all fields except change, which should be explicitly specified.
Example Usage:

An instance of InterpolParams can be created and used to control the bisection process within a line search method. For example:
```julia
interpolParams = InterpolParams(j=2, alphaj=50.0, change="decrease")
```
Ensure to correctly initialize and update all fields within these structures as per the requirements of your optimization algorithm.

"""
mutable struct InterpolParams
    j::Int
    alphaj::Float64
    alphaLo::Float64
    alphaHi::Float64
    alphatol::Float64
    alphatolBreached::Bool
    change::String

    function InterpolParams(;j=1, alphaj=100.0, alphaLo=0.0, alphaHi=100.0, alphatol=1e-10, alphatolBreached=false, change="noChange")
        new(j, alphaj, alphaLo, alphaHi, alphatol, alphatolBreached, change)
    end

end

"""
    CGStateType

A mutable struct designed to encapsulate the state and parameters pertinent to the Conjugate Gradient (CG) optimization method.

Fields:
- `k`: The global iteration counter across all optimization methods.
- `kCGD`: The iteration counter specific to the CG method.
- `gkm1`: The gradient vector at the previous CG iteration (`kCGD-1`).
- `gk`: The gradient vector at the current CG iteration (`kCGD`).
- `pkm1`: The search direction used in the previous CG iteration (`kCGD-1`).
- `pk`: The search direction for the current CG iteration (`kCGD`).
- `betakm1`: The β parameter value from the previous CG iteration used to compute `pkm1`.
- `betak`: The β parameter for the current CG iteration used to compute `pk`.
- `justRestarted`: A Boolean flag to indicate whether the CG method was restarted in the last CG iteration.

Constructor:
- The `CGStateType` constructor initializes the struct with the provided values or with default values if none are given. Default values are 1 for `k` and `kCGD`, zero vectors for `gkm1` and `gk`, zero vectors for `pkm1` and `pk`, and 0.0 for `betakm1` and `betak`. The `justRestarted` flag defaults to `false`.

Example Usage:
- To instantiate a `CGStateType` with default values:
    ```julia
    cgState = CGStateType()
    ```
- To set up with specific initial values for the CG method:
    ```julia
    cgState = CGStateType(k=3, kCGD=2, gkm1=[-1.0, -1.0], gk=[-0.5, -0.5], betakm1=0.5, betak=0.3, justRestarted=true)
    ```

This struct is utilized in each step of the CG method to maintain the necessary data for computing new search directions and for deciding when to restart the algorithm based on the `justRestarted` flag.
"""
mutable struct CGStateType
    k::Int
    kCGD::Int
    gkm1::Vector{Float64}
    gk::Vector{Float64}
    pkm1::Vector{Float64}
    pk::Vector{Float64}
    betakm1::Float64
    betak::Float64
    justRestarted::Bool

    function CGStateType(;
        k=1, kCGD=1,
        gkm1=Float64[], gk=Float64[],
        pkm1=Float64[], pk=Float64[], 
        betakm1=0.0, betak=0.0,
        justRestarted=false)

        new(k, kCGD, gkm1, gk, pkm1, pk, betakm1, betak, justRestarted)

    end

end

"""
    QNStateType

A mutable struct that holds the state of a Quasi-Newton (QN) optimization algorithm at each iteration.

Fields:
- `k`: The current iteration number.
- `xkm1`: The vector representing the solution at the previous iteration (`k-1`).
- `xk`: The vector representing the current solution at iteration `k`.
- `fkm1`: The value of the objective function at the previous iteration (`k-1`).
- `fk`: The value of the objective function at the current iteration (`k`).
- `gkm1`: The gradient vector of the objective function at the previous iteration (`k-1`).
- `gk`: The gradient vector of the objective function at the current iteration (`k`).
- `pkm1`: The search direction used in the previous iteration (`k-1`).
- `pk`: The search direction used in the current iteration (`k`).
- `Hkm1`: The approximation to the Hessian matrix or its inverse at the previous iteration (`k-1`).
- `Hk`: The approximation to the Hessian matrix or its inverse at the current iteration (`k`).

Constructor:
- The constructor for `QNStateType` can be called with keyword arguments for each field. If a field is not specified, it defaults to initial values such as `1` for `k`, `0.0` for scalar fields, empty arrays for vector fields, and an undefined `0x0` matrix for matrix fields.

Example Usage:
- To create a new instance of `QNStateType` with default initial values, you can simply call:
    ```julia
    qnState = QNStateType()
    ```
- To initialize with custom values for the current iteration, you could use:
    ```julia
    qnState = QNStateType(k=2, xkm1=[1.0, 1.0], xk=[1.5, 1.5], fkm1=10.0, fk=5.0)
    ```

This struct is instrumental in carrying forward information from one iteration to the next in a QN optimization routine, allowing for efficient updates of the solution vector and the approximation of the Hessian matrix.
"""
mutable struct QNStateType
    k::Int
    xkm1::Vector{Float64}
    xk::Vector{Float64}
    fkm1::Float64
    fk::Float64
    gkm1::Vector{Float64}
    gk::Vector{Float64}
    pkm1::Vector{Float64}
    pk::Vector{Float64}
    Hkm1::Matrix{Float64}
    Hk::Matrix{Float64}

    function QNStateType(;
        k=1, 
        xkm1=Float64[], xk=Float64[],
        fkm1=0.0, fk=0.0,
        gkm1=Float64[], gk=Float64[],
        pkm1=Float64[], pk=Float64[], 
        Hkm1=Matrix{Float64}(undef, 0, 0), Hk=Matrix{Float64}(undef, 0, 0))

        new(k, xkm1, xk, fkm1, fk, gkm1, gk, pkm1, pk, Hkm1, Hk)
    end

end

function TRparamsType(;
    Delta = 50.0,
    Delta_min = 1e-3,
    Delta_max = 100.0,
    etta_1 = 0.01,
    etta_2 = 0.25,
    etta_3 = 0.75,
    delta_1 = 0.25,
    delta_2 = 2.0,
    updateRadius = "uninitialized",
    accept = false
    )

    TRparams = Dict(
        :Delta => Delta,
        :Delta_min => Delta_min,
        :Delta_max => Delta_max,
        :etta_1 => etta_1,
        :etta_2 => etta_2,
        :etta_3 => etta_3,
        :delta_1 => delta_1,
        :delta_2 => delta_2,
        :updateRadius => updateRadius,
        :accept => accept
    )

    return TRparams
end

function SR1paramsType(;
    xkm1 = Vector{Float64}(undef, 0),
    gkm1 = Vector{Float64}(undef, 0),
    Bkm1 = Matrix{Float64}(undef, 0, 0)
    )

    SR1params = Dict(
        :xkm1 => xkm1,
        :gkm1 => gkm1,
        :Bkm1 => Bkm1
    )

    return SR1params
end

# SR1params = SR1paramsType()
# TRparams = TRparamsType()
# solState = SolStateType()
# solverState = SolverStateType()
# interpolParams = InterpolParams(alphatol=33)
# CGState = CGStateType()
# QNState = QNStateType()