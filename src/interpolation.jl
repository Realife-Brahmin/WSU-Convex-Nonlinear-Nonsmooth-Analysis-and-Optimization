using Parameters

include("types.jl")

"""
    bisection(interpolParams::InterpolParams)

Perform a bisection step to find a new trial step size `alphaj` within a specified interval.

# Arguments
- `interpolParams::InterpolParams`: A structure containing the parameters required for interpolation.

# Behavior
The function updates the bounds of the interval based on whether the step size needs to be increased or decreased.
The bisection process continues by incrementing the iteration number `j` and updating `alphaj`.

# Returns
- `interpolParams`: Updated interpolation parameters with the new trial step size and interval bounds.

# Examples
include("types.jl")

interpolParams = InterpolParams(j = 0, alphaj = 1.0, alphaLo = 0.0, alphaHi = 100.0, alphatol = 1e-8, change = "decrease")
updatedParams = bisection(interpolParams)

The `bisection` function is part of a larger optimization routine and expects that the calling function will handle the error conditions and the consequences of the `alphatolBreached` flag.
"""
function bisection(interpolParams::InterpolParams)
    @unpack j, alphaj, alphaLo, alphaHi, alphatol, change = interpolParams
    alphatolBreached = false

    if change == "increase"
        alphaLo = alphaj
    elseif change == "decrease"
        alphaHi = alphaj
    elseif change == "noChange"
        @error "Why is bisection even being performed if there's no change?"
    else
        @error "Bad condition"
    end

    j += 1
    alphaj = (alphaLo+alphaHi)/2
    if alphaHi - alphaLo < alphatol
        alphatolBreached = true
    end
    
    @pack! interpolParams = j, alphaj, alphaLo, alphaHi, alphatolBreached

    return interpolParams
end