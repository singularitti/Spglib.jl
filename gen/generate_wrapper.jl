using Clang

HEADER_BASE = "./deps/usr/include" # using system-wide installation of SCIP
all_headers = readdir(joinpath(HEADER_BASE, "spglib"))
context = Clang.init(
    headers = [joinpath(HEADER_BASE, "spglib", h) for h in all_headers],
    common_file = "spglib_common.jl",
    output_dir = "./src/Wrapper/",
    output_file = "spglib_api.jl",
    clang_args = ["-I", HEADER_BASE],
    header_wrapped = (header, cursorname) -> header == cursorname,
)
Clang.run(context)
