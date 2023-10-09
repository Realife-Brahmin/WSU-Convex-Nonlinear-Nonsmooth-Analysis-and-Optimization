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



