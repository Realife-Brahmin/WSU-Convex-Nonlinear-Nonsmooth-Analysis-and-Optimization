import Base.eval

function safe_import(mod_name::Symbol)
    eval(:(import $mod_name))
end

function safe_using(mod_name::Symbol)
    eval(:(using $mod_name))
end

# safe_using(:Objective)


