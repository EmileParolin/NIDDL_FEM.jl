# To build the documentation do
# $> julia --project=@. make.jl

push!(LOAD_PATH,"../src/")
using Documenter, NIDDL_FEM

makedocs(
    modules = [NIDDL_FEM],
    clean = false,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename="NIDDL_FEM.jl",
    authors="Emile Parolin",
)

deploydocs(
    repo = "github.com/EmileParolin/NIDDL_FEM.jl.git",
)
