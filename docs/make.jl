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
             "Simulation" => "sim.md"],
    format = Documenter.HTML(assets = ["assets/flux.css"])
)

 # deploydocs(
 #     repo = "github.com/SimonEnsemble/Controlz.jl.git",
 #     # This is a link to the main repo and the master branch
 #     # target = "build",
 #     julia = "1.3",
 #     osname = "linux",
 #     deps = Deps.pip("mkdocs", "mkdocs-material", "pymdown-extensions") # These are dependencies for the site, not the package
 # )
