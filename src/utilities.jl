# import Base.eval

function safe_import(mod_name::Symbol)
    eval(:(import $mod_name))
end

function safe_using(mod_name::Symbol)
    eval(:(using $mod_name))
end

# Extract fields from the 'res' structure and assign them to local variables.
# The purpose of this macro is to avoid repetitive assignments.
# Usage: @unpack_vars res status fvals αvals backtrackVals xvals
macro unpack_vars(struct_var, fields...)
    quote
        # For each specified field, assign the field's value from the structure to a local variable.
        $(Expr(:block, [:(local $(field) = $(struct_var).$(field)) for field in fields]...))
    end
end

# safe_using(:Objective)


