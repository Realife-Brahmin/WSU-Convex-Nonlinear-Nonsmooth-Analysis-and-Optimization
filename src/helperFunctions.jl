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

# Sample data
fvals = [i for i in 1:10]
gmagvals = [i^2 for i in 1:10]
αvals = [i^0.5 for i in 1:10]
backtrackVals = [10-i for i in 1:10]

xvals = hcat(1:10, 11:20)
gvals = hcat(21:30, 31:40)

data_1D = (fvals, gmagvals, αvals, backtrackVals)
data_2D = (xvals, gvals)

# Trimming function for 1D arrays
function trim_arrays(itr, data_tuple::Tuple...)
    return tuple([arr[1:itr-1] for arr in data_tuple]...)
end

# Trimming function for 2D arrays
function trim_arrays_2D(itr, data_tuple::Tuple...)
    return tuple([arr[:, 1:itr-1] for arr in data_tuple]...)
end

# Example usage with itr = 5
trimmed_data_1D = trim_arrays(5, data_1D)
trimmed_data_2D = trim_arrays_2D(5, data_2D)

println("Trimmed data_1D: ", trimmed_data_1D)
println("Trimmed data_2D: ", trimmed_data_2D)



