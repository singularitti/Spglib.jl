using Clang

const LIB_INCLUDE = "deps/usr/include"
const LIB_HEADERS = [joinpath(LIB_INCLUDE, "spglib", h) for h in readdir(joinpath(LIB_INCLUDE, "spglib"))]
context = Clang.init(
    headers = LIB_HEADERS,
    common_file = "spglib_common.jl",
    output_dir = "src/Wrapper/",
    output_file = "spglib_api.jl",
    clang_args = ["-I", LIB_INCLUDE],
    header_wrapped = (root, current) -> root == current,
    clang_diagnostics = true,
)
Clang.run(context)
