export get_error_code, get_error_message

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

get_error_code() = ccall((:spg_get_error_code, libsymspg), SpglibReturnCode, ())

get_error_message(spglib_error::SpglibReturnCode) = unsafe_string(
    ccall((:spg_get_error_message, libsymspg), Cstring, (SpglibReturnCode,), spglib_error),
)
