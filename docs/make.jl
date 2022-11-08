using AFMsimulations
using Documenter

DocMeta.setdocmeta!(AFMsimulations, :DocTestSetup, :(using AFMsimulations); recursive=true)

makedocs(;
    modules=[AFMsimulations],
    authors="Hannes Wallner",
    repo="https://github.com/HaneWall/AFMsimulations.jl/blob/{commit}{path}#{line}",
    sitename="AFMsimulations.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://HaneWall.github.io/AFMsimulations.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/HaneWall/AFMsimulations.jl",
    devbranch="main",
)
