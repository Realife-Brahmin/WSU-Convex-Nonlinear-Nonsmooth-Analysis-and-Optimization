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

function empty_FuncParam()::FuncParam
    return (params = Float64[], data = Matrix{Float64}(undef, 0, 0))
end

macro varname(var)
    return string(var)
end

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



