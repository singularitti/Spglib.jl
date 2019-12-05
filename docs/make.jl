using Documenter, Spglib

makedocs(;
    modules=[Spglib],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/singularitti/Spglib.jl/blob/{commit}{path}#L{line}",
    sitename="Spglib.jl",
    authors="Qi Zhang <singularitti@outlook.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/singularitti/Spglib.jl",
)
