using Parameters

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

## under testing
"""
    CGStateType

A mutable struct that encapsulates the state of the Conjugate Gradient (CG) method during optimization iterations.

Fields:
- `k`: The current iteration number in the CG method.
- `xkm1`: The solution vector at the previous iteration (k-1).
- `xk`: The current solution vector at iteration k.
- `fkm1`: The value of the objective function at the previous iteration (k-1).
- `fk`: The current value of the objective function at iteration k.
- `gkm1`: The gradient vector of the objective function at the previous iteration (k-1).
- `gk`: The current gradient vector of the objective function at iteration k.
- `gmagkm1`: The magnitude of the gradient at the previous iteration (k-1).
- `gmagk`: The current magnitude of the gradient at iteration k.
- `pkm1`: The search direction used in the previous iteration (k-1).
- `pk`: The current search direction at iteration k.
- `betakm1`: The beta coefficient used to update the search direction in the previous iteration (k-1).
- `betak`: The current beta coefficient used to update the search direction at iteration k.
- `justRestarted`: A Boolean flag indicating if the CG method was just restarted.

Constructor:
The constructor can be invoked with keyword arguments for all the fields. Default values are used when specific values are not provided by the caller.

Example Usage:
Instantiate the CGStateType with the desired initial values or use the defaults for starting a new CG optimization routine.
```julia
cgState = CGStateType(k=2, xkm1=[0.0, 0.0], xk=[1.0, 2.0], fkm1=5.0, fk=3.0)
```
This instance `cgState` now contains all the necessary state information for the second iteration of a CG optimization process. The fields `xkm1` and `xk` are initialized to vectors representing positions in the parameter space, `fkm1` and `fk` are initialized to the respective objective function values, and all other vector and scalar fields are set according to the provided values. The justRestarted flag is set to `false`, indicating that the CG method is continuing from the previous state.
"""
mutable struct CGStateType
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
    betakm1::Float64
    betak::Float64
    justRestarted::Bool

    function CGStateType(;
        k=1, 
        xkm1=Float64[], xk=Float64[],
        fkm1=0.0, fk=0.0,
        gkm1=Float64[], gk=Float64[],
        gmagkm1=0.0, gmagk=0.0, 
        pkm1=Float64[], pk=Float64[], 
        betakm1=0.0, betak=0.0,
        justRestarted=false)

        new(k, xkm1, xk, fkm1, fk, gkm1, gk, gmagkm1, gmagk, pkm1, pk, betakm1, betak, justRestarted)
    end

end

# under testing
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

# solState = SolStateType()
# solverState = SolverStateType()
# interpolParams = InterpolParams(alphatol=33)
# CGState = CGStateType()
# QNState = QNStateType()