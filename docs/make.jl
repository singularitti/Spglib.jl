using Spglib
using Documenter

DocMeta.setdocmeta!(Spglib, :DocTestSetup, :(using Spglib); recursive=true)

makedocs(;
    modules=[Spglib],
    authors="Qi Zhang <singularitti@outlook.com>",
    repo="https://github.com/singularitti/Spglib.jl/blob/{commit}{path}#{line}",
    sitename="Spglib.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://singularitti.github.io/Spglib.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Installation" => "installation.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/singularitti/Spglib.jl",
)
