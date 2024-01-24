export get_error_code, get_error_message

# See https://www.mcobject.com/docs/Content/Programming/C/Return_Codes.htm or
# https://www.gnu.org/software/libc/manual/html_node/Exit-Status.html
"Represent various return codes from the Spglib library."
@enum SpglibReturnCode::Int8 begin
    SUCCESS = 0
    SPACEGROUP_SEARCH_FAILED = 1
    CELL_STANDARDIZATION_FAILED = 2
    SYMMETRY_OPERATION_SEARCH_FAILED = 3
    ATOMS_TOO_CLOSE = 4
    POINTGROUP_NOT_FOUND = 5
    NIGGLI_FAILED = 6
    DELAUNAY_FAILED = 7
    ARRAY_SIZE_SHORTAGE = 8
    NONE = 9
end

struct SpglibError <: Exception
    msg::String
end

# The return value for the C function `spg_get_error_code` is a C enum called
# `SpglibError`. The C compiler guarantees that enums are internally stored as
# type `int`, which Julia represents as `Cint`. The size of `int` is
# "implementation dependent". Frequently 4 or 8 bytes, but this depends on
# architecture and compiler details: https://stackoverflow.com/a/366033/500314.
# Directly casting this Cint into an @enum or SumTypes.jl value would be memory
# unsafe (see https://github.com/singularitti/Spglib.jl/issues/183).
"""
    get_error_code()

Return an instance of the enumerated type `SpglibReturnCode`.
"""
get_error_code() = SpglibReturnCode(@ccall libsymspg.spg_get_error_code()::Cint)

"""
    get_error_message(code::SpglibReturnCode)

Get the corresponding error message for a given error code.
"""
get_error_message(code::SpglibReturnCode) =
    unsafe_string(@ccall libsymspg.spg_get_error_message(Cint(code)::Cint)::Cstring)

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
