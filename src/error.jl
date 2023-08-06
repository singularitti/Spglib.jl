export get_error_code, get_error_message, check_error

# See https://www.mcobject.com/docs/Content/Programming/C/Return_Codes.htm or
# https://www.gnu.org/software/libc/manual/html_node/Exit-Status.html
@enum SpglibReturnCode begin
    SUCCESS  # 0
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

get_error_code() = ccall((:spg_get_error_code, libsymspg), SpglibReturnCode, ())

get_error_message(code::SpglibReturnCode) = unsafe_string(
    ccall((:spg_get_error_message, libsymspg), Cstring, (SpglibReturnCode,), code)
)

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
