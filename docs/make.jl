using Spglib
using Documenter

DocMeta.setdocmeta!(Spglib, :DocTestSetup, :(using Spglib); recursive=true)

makedocs(;
    modules=[Spglib],
    authors="singularitti <singularitti@outlook.com> and contributors",
    repo="https://github.com/singularitti/Spglib.jl/blob/{commit}{path}#{line}",
    sitename="Spglib.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://singularitti.github.io/Spglib.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Installation Guide" => "man/installation.md",
            "Troubleshooting" => "man/troubleshooting.md",
        ],
        "Reference" => Any[
            "Public API" => "lib/public.md",
            # "Internals" => map(
            #     s -> "lib/internals/$(s)",
            #     sort(readdir(joinpath(@__DIR__, "src/lib/internals")))
            # ),
        ],
        "Developer Docs" => [
            "Contributing" => "developers/contributing.md",
            "Style Guide" => "developers/style-guide.md",
            "Design Principles" => "developers/design-principles.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/singularitti/Spglib.jl",
    devbranch="main",
)
