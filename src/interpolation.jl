using Parameters

include("types.jl")

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