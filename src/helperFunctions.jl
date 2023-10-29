"""
    myprintln(print_flag::Bool, message; log::Bool=true, log_path::String="./logging/logs.txt")

Print a message to the console if `print_flag` is `true`. Additionally, if `log` is `true`,
the message will also be written to a log file specified by `log_path`. By default, the log path
is set to `"./logging/logs.txt"`.
"""
function myprintln(print_flag::Bool, message; log::Bool=true, log_path::String="./logging/logs.txt")
    if print_flag
        println(message)
        if log
            open(log_path, "a") do f  # append mode
                println(f, message)
            end
        end
    end
end


"""
    empty_FuncParam()::FuncParam

Return an empty `FuncParam` object. The returned object has `params` set to an empty Float64 vector
and `data` set to an empty Float64 matrix.
"""
function empty_FuncParam()::FuncParam
    return (params = Float64[], data = Matrix{Float64}(undef, 0, 0))
end


"""
    @checkForNaN(varname)

A macro that checks a vector for any NaN values. If NaN values are found,
prints a message indicating the vector's name and the index where NaN is located.
"""
macro checkForNaN(varname)
    vec = esc(varname)
    quote
        for (idx, val) in enumerate($vec)
            if isnan(val)
                println("Vector ", string($(QuoteNode(varname))), " has NaN at index ", idx)
            end
        end
    end
end


"""
    myzeros(array)

Create a new array with the same size and type as the given `array`, but initialized with zeros.
"""
function myzeros(array)
    return zeros(eltype(array), size(array))
end

"""
    myfill(array, value)

Create a new array with the same size and type as the given `array`, but initialized with `value`.
The type of the returned array's elements is preserved from the original array.
"""
function myfill(array, value)
    T = eltype(array)
    return fill(T(value), size(array))
end


