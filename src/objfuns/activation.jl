function activation(x; 
    functionName::String = "sigmoid")

    if functionName == "sigmoid"
        y = 1.0/( 1.0 + exp(-x) )
    else
        @error "Other activation functions currently not supported"
    end

    return y
end

# activation(-1)
# activation(-0.5)
# activation(0.0)
# activation(1)