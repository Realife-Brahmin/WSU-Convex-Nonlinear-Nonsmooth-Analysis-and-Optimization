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

"""
    trim_array(nt::NamedTuple, n::Int) -> NamedTuple

Trim arrays within a NamedTuple to a predetermined maximum length or size.

# Arguments
- `nt::NamedTuple`: A NamedTuple containing arrays as values.
- `n::Int`: Maximum size for arrays. 1D arrays exceeding this length will be trimmed. 
For 2D arrays, both rows and columns will be trimmed if they exceed this size.

# Returns
- A new NamedTuple with trimmed arrays.

# Example
nt = (a = [1, 2, 3, 4, 5], b = [6, 7, 8, 9, 10], c = rand(10, 10), d = 11)
trimmed_nt = trim_array(nt, 3)
println(trimmed_nt)

# Notes
- This function is designed to handle 1D (vectors) and 2D arrays (matrices). 
Other data structures remain unchanged.
- It utilizes a dictionary as an intermediary data structure before converting back to a NamedTuple.
"""
function trim_array(nt::NamedTuple, n::Int)
    # Create a dictionary to hold the trimmed arrays
    trimmed_dict = Dict()

    # Iterate over each key in the NamedTuple
    for key in keys(nt)
        # Get the array using getproperty
        arr = getproperty(nt, key)
        
        # Check if it's a 1D array and trim it
        if isa(arr, AbstractVector) && length(arr) > n
            arr = arr[1:n]
        # Check if it's a 2D array and trim it
        elseif isa(arr, AbstractMatrix) && any(size(arr) .> n)
            arr = arr[1:min(end,n), 1:min(end,n)]
        end

        # Store the trimmed array in the dictionary
        trimmed_dict[key] = arr
    end
    
    # Convert the dictionary back to a NamedTuple and return
    return NamedTuple(trimmed_dict)
end





