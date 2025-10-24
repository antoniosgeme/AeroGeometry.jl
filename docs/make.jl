using Documenter
using AeroGeometry

makedocs(
    sitename = "AeroGeometry.jl",
    format = Documenter.HTML(),
    modules = [AeroGeometry],
    pages = [
        "Home" => "index.md",
    ]
)

deploydocs(
    repo = "github.com/antoniosgeme/AeroGeometry.jl.git",
    devbranch = "master"
)
