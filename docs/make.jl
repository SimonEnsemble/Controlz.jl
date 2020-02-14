using Pkg;
Pkg.activate(joinpath(@__DIR__, "..")); Pkg.instantiate()
Pkg.activate(); Pkg.instantiate()

pushfirst!(LOAD_PATH, joinpath(@__DIR__, ".."))
# ^ above from Flux.jl to get this to work

using Documenter
using Controlz

makedocs(
    root = joinpath(dirname(pathof(Controlz)), "..", "docs"),
    modules = [Controlz],
 #     source = "src",
    sitename = "Controlz.jl",
    clean = true,
    pages = ["Controlz" => "index.md",
             "Transfer Functions" => "tfs.md",
             "Simulation" => "sim.md",
             "Visualization" => "viz.md",
             "Control systems" => "controls.md",
             "FAQ" => "faq.md"]
 #     format = Documenter.HTML(assets = ["assets/flux.css"])
)

deploydocs(
    repo = "github.com/SimonEnsemble/Controlz.jl.git",
    # This is a link to the main repo and the master branch
    # target = "build",
    deps = Deps.pip("mkdocs", "mkdocs-bootswatch", "pymdown-extensions") # These are dependencies for the site, not the package
)
