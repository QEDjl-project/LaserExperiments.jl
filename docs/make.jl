using Pkg

# targeting the correct source code
# this assumes the make.jl script is located in LaserExperiments.jl/docs
project_path = Base.Filesystem.joinpath(Base.Filesystem.dirname(Base.source_path()), "..")
Pkg.develop(; path = project_path)

using Documenter
using LaserExperiments

DocMeta.setdocmeta!(LaserExperiments, :DocTestSetup, :(using LaserExperiments); recursive = true)

# some paths for links
readme_path = joinpath(project_path, "README.md")
index_path = joinpath(project_path, "docs/src/index.md")
license_url = "https://github.com/QEDjl-project/LaserExperiments.jl/blob/main/LICENSE"
code_of_conduct_url = "https://github.com/QEDjl-project/LaserExperiments.jl/blob/main/CODE_OF_CONDUCT.md"

# Copy README.md from the project base folder and use it as the start page
open(readme_path, "r") do readme_in
    readme_string = read(readme_in, String)

    # replace relative links in the README.md
    readme_string = replace(readme_string, "[MIT](LICENSE)" => "[MIT]($(license_url))")
    readme_string = replace(readme_string, "(CODE_OF_CONDUCT.md)" => "($(code_of_conduct_url))")
    readme_string = replace(readme_string, "[contributing guide directly on GitHub](docs/src/90-contributing.md)" => "contributing guide directly on GitHub")

    open(index_path, "w") do readme_out
        write(readme_out, readme_string)
    end
end

pages = [
    "Home" => "index.md",
    "Contributing" => "90-contributing.md",
    "Dev Docs" => "91-developer.md",
    "References" => "95-reference.md",
]

try
    # generate docs with Documenter.jl
    makedocs(;
        modules = [LaserExperiments],
        checkdocs = :exports,
        authors = "Uwe Hernandez Acosta",
        repo = Documenter.Remotes.GitHub("QEDjl-project", "LaserExperiments.jl"),
        sitename = "LaserExperiments.jl",
        format = Documenter.HTML(;
            prettyurls = get(ENV, "CI", "false") == "true",
            canonical = "https://qedjl-project.gitlab.io/LaserExperiments.jl",
            assets = String[],
            mathengine = Documenter.MathJax2(),
            collapselevel = 1,
        ),
        pages = pages,
    )
finally
    # doing some garbage collection
    @info "GarbageCollection: remove generated landing page"
    rm(index_path)
end

deploydocs(; repo = "github.com/QEDjl-project/LaserExperiments.jl.git", push_preview = false)
