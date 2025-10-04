using LaserExperiments
using Documenter

DocMeta.setdocmeta!(LaserExperiments, :DocTestSetup, :(using LaserExperiments); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
        file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [LaserExperiments],
    authors = "Uwe Hernandez Acosta <u.hernandez@hzdr.de>",
    repo = Documenter.Remotes.GitHub("QEDjl-project", "LaserExperiments.jl"),
    sitename = "LaserExperiments.jl",
    #format = Documenter.HTML(; canonical = "https://QEDjl-project.github.io/LaserExperiments.jl"),
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://qedjl-project.gitlab.io/LaserExperiments.jl",
        assets = String[],
        mathengine = Documenter.MathJax2(),
        collapselevel = 1,
        edit_link = :commit,
        # TODO: workaround
        # should be fixed: https://github.com/QEDjl-project/QEDbase.jl/issues/4
        size_threshold_ignore = ["index.md"],
    ),
    pages = ["index.md"; numbered_pages],
)

deploydocs(;
    devbranch = "main",
    repo = "github.com/QEDjl-project/LaserExperiments.jl",
    push_preview = true
)
