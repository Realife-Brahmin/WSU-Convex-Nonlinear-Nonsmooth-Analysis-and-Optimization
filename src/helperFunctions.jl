"""
    @packDict(vars_str)

Create a dictionary from a string of variable names formatted as `"{var1, var2, var3, ...}"`.

# Arguments
- `vars_str::String`: A single string argument formatted with variable names enclosed in curly braces and separated by commas.

# Returns
- `Dict{Symbol, Any}`: Returns a dictionary where each key is a symbol corresponding to a variable name, and each value is the current value of that variable at the time the macro is invoked.

# Usage
The macro automates the creation of a dictionary by parsing a string containing variable names, which are separated by commas and enclosed in curly braces. It constructs a dictionary where keys are the variable names as symbols and values are their corresponding current values in the scope where the macro is invoked.

# Notes
- **Security**: Use caution with `eval()` as it will execute code with all available permissions, which can be a security risk if variable contents are not controlled.
- **Scope**: Variables must be defined in the scope where the macro is called; otherwise, `eval()` will not be able to find and return their values.
- **Performance**: Since `eval()` is used, the performance may be impacted if the macro is used in performance-critical parts of your application.
- **Mutation Warning**: If a variable is modified (mutated) after being initially defined, the `@packDict` macro will still pack the original value from when it was first introduced into the scope, unless the variable is reassigned. Reassignment must occur in the same or a broader scope than where `@packDict` is invoked.

# Examples
```julia
# Define some variables
a = 10
b = [1, 2, 3]
c = "Hello, world!"

# Create a dictionary from these variables using the macro
myDict = @packDict "{a, b, c}"

# Output the created dictionary
println("Dictionary contents: ", myDict)

# Expected Output
Dictionary contents: Dict(:a => 10, :b => [1, 2, 3], :c => "Hello, world!")
```
"""
macro packDict(vars_str)
    # Remove the curly braces and split the string into symbols
    vars = strip(vars_str[2:end-1])  # Strip out the enclosing braces
    vars_list = split(vars, ",")     # Split the string into components based on commas
    symbol_list = [Symbol(strip(v)) for v in vars_list]  # Convert strings to symbols and strip whitespace

    # Generate the dictionary packing expression
    esc(quote
        dict_expr = Dict{Symbol, Any}()
        for var in $(symbol_list)  # Ensure symbol_list is treated as a literal array of symbols
            dict_expr[var] = eval(var)
        end
        dict_expr
    end)
end


"""
    initializeLogFile(logPath::String)

Initializes the logging environment by ensuring that the default log file is cleared at the start of the program. If the log file already exists, it is removed to start fresh, facilitating clean and relevant logging for the current execution.

# Arguments
- `logPath`: The file path to the log file. If the file exists, it will be removed.

# Usage Example
```julia
# Specify the default log file path
logDefault = "./logging/logs.txt"

# Initialize the log file, clearing any existing content
initializeLogFile(logDefault)
"""
function initializeLogFile(logPath::String = "./logging/logs.txt")
    if isfile(logPath)
        rm(logPath)
        myprintln(true, "Existing log file at '$logPath' removed.", color=:yellow)
    else
        myprintln(true, "No existing log file found at '$logPath'. Starting fresh.", color=:yellow)
    end
end

"""
    myprintln(print_flag::Bool, message; color=:normal, log::Bool=true, log_path::String="./logging/logs.txt")

Print a styled message to the console if `print_flag` is `true`. If `color` is specified, the message will
be printed in the specified color. Additionally, if `log` is `true`, the message will also be written to a 
log file specified by `log_path` without color styling. By default, the log path is set to `"./logging/logs.txt"`.
"""
function myprintln(print_flag::Bool, message;
    color=:normal,
    log::Bool=true,
    log_path::String="./logging/logs.txt"
)

    if print_flag
        printstyled(message; color=color)  # Apply color styling
        println()  # Move to the next line
        if log
            open(log_path, "a") do f  # append mode
                println(f, message)  # Log the message without color information
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


# function trim_array(nt::NamedTuple, itrMax::Int;
#     keys_to_remove::Vector{Symbol}=[:xopt, :fopt])
#     # Create a dictionary to hold the trimmed arrays
#     trimmed_dict = Dict()

#     # Convert the keys of the named tuple to a set
#     keys_set = Set(keys(nt))

#     # Define the keys you want to remove as another set
#     keys_to_remove = Set(keys_to_remove)    

#     # Use set difference to remove the specified keys
#     @show filtered_keys_set = setdiff(keys_set, keys_to_remove)

#     # If you need to work with a vector of keys afterwards
#     @show filtered_keys_vector = collect(filtered_keys_set)

#     # Iterate over each key in the NamedTuple
#     for key in filtered_keys_vector
#         # Get the array using getproperty
#         arr = getproperty(nt, key)
        
#         # @show key
#         # Check if it's a 1D array and trim it
#         if eltype(arr) ∈ (Int64, Float64) && isa(arr, AbstractVector)
#             arr = arr[1:itrMax]
#         # Check if it's a 2D array and trim it
#         elseif eltype(arr) ∈ (Int64, Float64) && isa(arr, AbstractMatrix)
#             arr = arr[:, 1:itrMax]
#         end

#         # Store the trimmed array in the dictionary
#         trimmed_dict[key] = arr
#     end
    
#     # Convert the dictionary back to a NamedTuple and return

#     @show trimmed_dict[:xopt], trimmed_dict[:fopt]
#     return NamedTuple(trimmed_dict)
# end

function trim_array(nt::NamedTuple, itrMax::Int;
    keys_to_not_trim::Vector{Symbol}=[:xopt, :fopt])
    # Create a dictionary to hold the possibly modified values
    trimmed_dict = Dict()

    # Convert the keys of the named tuple to a set
    keys_set = Set(keys(nt))

    # Define the keys you want to potentially trim as another set
    keys_to_trim_set = setdiff(keys_set, keys_to_not_trim)

    # Iterate over each key in the NamedTuple
    for key in keys(nt)
        # Get the current value using getproperty
        value = getproperty(nt, key)

        # Decide if the key's value needs trimming
        if key in keys_to_trim_set
            # Check if it's a 1D array and trim it
            if isa(value, AbstractVector) && eltype(value) in (Int64, Float64)
                value = value[1:itrMax]
                # Check if it's a 2D array and trim it
            elseif isa(value, AbstractMatrix) && eltype(value) in (Int64, Float64)
                value = value[:, 1:itrMax]
            end
        end

        # Store the original or trimmed value in the dictionary
        trimmed_dict[key] = value
    end
        # @show trimmed_dict[:xopt], trimmed_dict[:fopt]

    # Convert the dictionary back to a NamedTuple and return
    return NamedTuple(trimmed_dict)
end


"""
    extrapolate(xopt::Vector{Float64}, factor::Int) -> Vector{Float64}

Generate a new vector that extrapolates a set of monotonically increasing values `xopt` by inserting `factor - 1` evenly spaced points between each consecutive pair of values. The new vector will have `factor * length(xopt)` elements.

# Arguments
- `xopt::Vector{Float64}`: The original vector of monotonically increasing values, not including `0.0` but including the end point `1.0`.
- `factor::Int`: The number of points to be included between each pair of the original vector, including one of the bounds.

# Returns
- `Vector{Float64}`: A new vector containing the original and the inserted points.

# Example
```julia
xopt = [0.1, 0.5, 0.7]
factor = 3  # This will insert 2 evenly spaced values between each pair in xopt
x0 = extrapolate(xopt, factor)
# Output: [0.1, 0.233..., 0.366..., 0.5, 0.6, 0.7, 0.8, 0.9]
```
The function ensures that the new points are distributed evenly between the existing points from xopt, with the spacing determined by the factor. The original points in xopt are included in the output, and the interpolation considers the natural boundary at 1.0.

## Note
This function assumes that xopt does not include the boundary point 0.0 and that the last point is less than 1.0. The interpolation will add points up to the boundary at 1.0.
"""
function extrapolate(xopt::Vector{Float64}, factor::Int)
    n = length(xopt)
    # The new array will have the original points plus (factor - 1) points between each
    new_length = n * factor
    new_xopt = Vector{Float64}(undef, new_length)

    # Fill the first element with xopt[1]
    new_xopt[1] = xopt[1]
    
    # Fill in the rest of the points
    for i in 1:n-1
        # Calculate the interpolated values
        for j in 1:factor-1
            new_xopt[(i-1)*factor + j + 1] = ((factor - j) * xopt[i] + j * xopt[i + 1]) / factor
        end
        # Place the original value
        new_xopt[i*factor + 1] = xopt[i + 1]
    end
    
    # Handle the last segment between xopt[end] and 1.0
    for j in 1:factor-1
        new_xopt[(n-1)*factor + j + 1] = ((factor - j) * xopt[end] + j) / factor
    end
    
    return new_xopt
end

import Interpolations
"""
    interpolate_velocity(v::Array{Float64,2}, xxm::Vector{Float64}, yym::Vector{Float64}) -> Vector{Float64}

Performs bilinear interpolation on a 2D grid of velocity data.

# Arguments
- `v::Array{Float64,2}`: A 2D array representing velocity values on a grid. The array can be of any rectangular shape.
- `xxm::Vector{Float64}`: A vector containing the x-coordinates at which to interpolate the velocity. The values should correspond to the indices of the grid and can be floating-point numbers.
- `yym::Vector{Float64}`: A vector containing the y-coordinates at which to interpolate the velocity, similar to `xxm`.

# Returns
- `Vector{Float64}`: A vector of interpolated velocity values corresponding to each pair of coordinates in `xxm` and `yym`.

# Notes
- The function uses bilinear interpolation, which is suitable for continuous and smoothly varying velocity fields.
- The length of `xxm` and `yym` must be the same, as they represent pairs of coordinates in the velocity grid.
- This function scales the interpolation object to match the index coordinates of the original matrix, allowing for direct use of the original floating-point coordinates.

# Examples
To use this function, ensure that `v` is a 2D array of velocity values, and `xxm` and `yym` are vectors containing the coordinates for interpolation. The function will return a vector of interpolated values.

```julia
v = [sin(i + j) for i in 1:256, j in 1:300]
xxm = rand(1:256, 1000) 
yym = rand(1:300, 1000)
vel = interpolate_velocity(v, xxm, yym)
```
"""
function interpolate_velocity(v, xxm, yym)
    # Check if the length of xxm and yym are the same
    if length(xxm) != length(yym)
        error("Length of xxm and yym must be the same")
    end

    # Create an interpolation object for bilinear interpolation
    interp_v = Interpolations.interpolate(v, Interpolations.BSpline(Interpolations.Linear()), Interpolations.OnGrid())

    # Scale the interpolation object to match the index coordinates of the original matrix
    scaled_interp_v = Interpolations.scale(interp_v, 1:size(v, 1), 1:size(v, 2))

    # Perform interpolation
    [scaled_interp_v[x, y] for (x, y) in zip(xxm, yym)]
end;




