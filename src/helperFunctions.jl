function myprintln(verbose::Bool, message::Any)
    if verbose
        println(message)
    end
end