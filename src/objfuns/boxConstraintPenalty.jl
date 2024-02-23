function boxConstraintPenalty(x::Vector, box::Dict)

    @unpack indices, lbs, ubs, barrier = box

    P = 0.0
    for constraintNum ∈ eachindex(indices)
        i = indices[constraintNum]
        P += barrier * (min(0, x[i] - lbs[constraintNum])^2)
        P += barrier * (min(0, ubs[constraintNum] - x[i])^2)
    end
    return P
end

function boxConstraintGradient(x::Vector, box::Dict)

    @unpack indices, lbs, ubs, barrier = box

    gP = 0.0

    for constraintNum ∈ eachindex(indices)
        i = indices[constraintNum]
        gP += -2 * barrier * min(0, x[i] - lbs[constraintNum])
        gP += -2 * barrier * min(0, ubs[constraintNum] - x[i])
    end

    return gP
end