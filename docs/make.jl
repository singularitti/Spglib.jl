using Documenter, SpgLib

makedocs(;
    modules=[SpgLib],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/singularitti/SpgLib.jl/blob/{commit}{path}#L{line}",
    sitename="SpgLib.jl",
    authors="Qi Zhang <singularitti@outlook.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/singularitti/SpgLib.jl",
)
