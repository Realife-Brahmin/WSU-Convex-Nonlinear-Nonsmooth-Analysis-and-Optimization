# import Base.eval

function safe_import(mod_name::Symbol)
    eval(:(import $mod_name))
end

function safe_using(mod_name::Symbol)
    eval(:(using $mod_name))
end

macro unpack_vars(namedtuple_var)
    exprs = []
    fields = fieldnames(typeof(namedtuple_var))
    for field in fields
        push!(exprs, :( $(Symbol(field)) = $namedtuple_var.$(Symbol(field)) ))
    end
    return esc(quote
        $(exprs...)
    end)
end

# a = "aa2"
# b = "baba2"
# res1 = (a=a, b=b)
# @unpack_vars res1
# a
# b
# res1.a
# res1.b

# safe_using(:Objective)


