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
             "FAQ" => "faq.md"],
    format = Documenter.HTML(assets = ["assets/flux.css"])
)

deploydocs(
    repo = "github.com/SimonEnsemble/Controlz.jl.git",
    # This is a link to the main repo and the master branch
    # target = "build",
    # versions sugggested by https://discourse.julialang.org/t/mkdocs-material-in-documenter/13764
    deps = Deps.pip("mkdocs==0.17.5", "mkdocs-material==2.9.4", "pymdown-extensions") # These are dependencies for the site, not the package
)
