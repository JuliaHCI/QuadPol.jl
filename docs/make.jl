using QuadPol
using Documenter

DocMeta.setdocmeta!(QuadPol, :DocTestSetup, :(using QuadPol); recursive=true)

makedocs(;
    modules=[QuadPol],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo="https://github.com/JuliaHCI/QuadPol.jl/blob/{commit}{path}#{line}",
    sitename="QuadPol.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaHCI.github.io/QuadPol.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaHCI/QuadPol.jl",
    devbranch="main",
)
