using Muspel
using Documenter

DocMeta.setdocmeta!(Muspel, :DocTestSetup, :(using Muspel); recursive=true)

makedocs(;
    modules=[Muspel],
    authors="Tiago M. D. Pereira",
    repo="https://github.com/tiagopereira/Muspel.jl/blob/{commit}{path}#{line}",
    sitename="Muspel.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tiagopereira.github.io/Muspel.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tiagopereira/Muspel.jl",
    devbranch="main",
)
