function round_to_n_digits(x::Vector{Float64}; digits=9)
    return round.(x, digits=digits)
    # return trunc.(x, digits=digits)
end

function simple_derivative(f, x::Vector{Float64}; δ::Float64=round_to_n_digits(sqrt(1e-9)))
    δ = round_to_n_digits(δ)
    xcap = round_to_n_digits(normalize(x))
    xk = round_to_n_digits(x)
    xkp1 = round_to_n_digits(xk + δ*xcap)
    fk = round_to_n_digits(f(xk))
    fkp1 = round_to_n_digits(f(xkp1))

    df = round_to_n_digits( (fkp1 - fk)/δ )
    return round_to_n_digits(df)
end

function inClassExample()
    ff = x -> sin.(x)
    gg = x -> cos.(x)

    δVals = round_to_n_digits.(vcat(10.0 .^ (-9:-5), sqrt(1e-9), 10.0 .^ (-4:0)))

    xk = [1.0]
    gk = round_to_n_digits(gg([1.0]))
    digits = 9
    myprintln(true, "*"^50)
    myprintln(true, "Here, the precision has been set to $(digits) digits.")
    keepRunningIterations = true
    errorMin = 100.0

    for δ = δVals
        myprintln(true, "*"^33)
        myprintln(true, "For step δ = $(δ)")
        ∇fk = simple_derivative(ff, xk, δ=δ)
        myprintln(true, "∇f = $(∇fk)") 
        global gk
        # println(gk)
        local err = maximum(abs.(gk-∇fk))
        myprintln(true, "Error = $(err)")
        global errorMin
        if err < errorMin
            errorMin = err
            global best_step = δ
        end
    end

    myprintln(true, "*"^45)
    myprintln(true, "Least error happened at δ=$(best_step) with error = $(errorMin)")
end