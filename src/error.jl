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
function get_error_code()
    # The return value for the C function `spg_get_error_code` is a C enum
    # called SpglibError. The C compiler guarantees that enums are internally
    # stored as type `int`; its size is "implementation dependent". Typically 4
    # or 8 bytes, but this could depend on the architecture and compiler details
    # https://stackoverflow.com/a/366033/500314. Julia provides `CInt` to
    # guarantee an integer that has the same size as the C `int` value. Manually
    # map this enum tag to variants of the `@enum SpglibReturnCode`. The latter
    # is defined as a Union using the Julia SumTypes library, and its
    # representation details are also implementation dependent (e.g., the "tag"
    # of this Union need not have the same size as a `CInt`). Note that a
    # mistake here could cause memory corruption, as discussed in
    # https://github.com/singularitti/Spglib.jl/issues/183.
    code = @ccall libsymspg.spg_get_error_code()::CInt
    if code == 0
        return SUCCESS
    elseif code == 1
        return SPACEGROUP_SEARCH_FAILED
    elseif code == 2
        return CELL_STANDARDIZATION_FAILED
    elseif code == 3
        return SYMMETRY_OPERATION_SEARCH_FAILED
    elseif code == 4
        return ATOMS_TOO_CLOSE
    elseif code == 5
        return POINTGROUP_NOT_FOUND
    elseif code == 6
        return NIGGLI_FAILED
    elseif code == 7
        return DELAUNAY_FAILED
    elseif code == 8
        return ARRAY_SIZE_SHORTAGE
    elseif code == 9
        return NONE
    end
end

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
