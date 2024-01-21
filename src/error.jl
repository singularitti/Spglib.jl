using SumTypes: @sum_type

export get_error_code, get_error_message

# See https://www.mcobject.com/docs/Content/Programming/C/Return_Codes.htm or
# https://www.gnu.org/software/libc/manual/html_node/Exit-Status.html
"Represent various return codes from the Spglib library."
@sum_type SpglibReturnCode begin
    SUCCESS
    SPACEGROUP_SEARCH_FAILED
    CELL_STANDARDIZATION_FAILED
    SYMMETRY_OPERATION_SEARCH_FAILED
    ATOMS_TOO_CLOSE
    POINTGROUP_NOT_FOUND
    NIGGLI_FAILED
    DELAUNAY_FAILED
    ARRAY_SIZE_SHORTAGE
    NONE
end

struct SpglibError <: Exception
    msg::String
end

"""
    get_error_code()

Return an instance of the enumerated type `SpglibReturnCode`.
"""
get_error_code() = @ccall libsymspg.spg_get_error_code()::SpglibReturnCode

"""
    get_error_message(code::SpglibReturnCode)

Get the corresponding error message for a given error code.
"""
get_error_message(code::SpglibReturnCode) =
    unsafe_string(@ccall libsymspg.spg_get_error_message(code::SpglibReturnCode)::Cstring)

"""
    check_error()

Check if an error has occurred and throw an exception if so.
"""
function check_error()
    code = get_error_code()
    if code == SUCCESS
        return nothing
    else
        return throw(SpglibError(get_error_message(code)))
    end
end

# See https://github.com/JuliaLang/julia/blob/3903fa5/base/missing.jl#L18-L19
Base.showerror(io::IO, ex::SpglibError) = print(io, "SpglibError: ", ex.msg, '!')
