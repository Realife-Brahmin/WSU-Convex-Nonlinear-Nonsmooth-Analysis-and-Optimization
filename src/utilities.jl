function replace_field(nt::NamedTuple, field::Symbol, value)
    fields = fieldnames(typeof(nt))
    values = map(fields) do f
        if f == field
            value
        else
            getproperty(nt, f)
        end
    end
    return NamedTuple{fields}(values)
end


