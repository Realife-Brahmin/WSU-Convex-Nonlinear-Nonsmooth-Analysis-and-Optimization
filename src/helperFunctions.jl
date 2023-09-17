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

function initialize_logging(;log_path::String="./logging/logs.txt", overwrite::Bool=false)
    if overwrite
        # Clear the old log by opening it in write mode and immediately closing it.
        open(log_path, "w") do f end
    end
end


