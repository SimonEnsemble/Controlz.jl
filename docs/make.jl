using Documenter
using Controlz
using Polynomials

# using Controlz before all docstrings
DocMeta.setdocmeta!(Controlz, :DocTestSetup, :(using Controlz); recursive=true)

makedocs(
    root = joinpath(dirname(pathof(Controlz)), "..", "docs"),
    modules = [Controlz],
 #     source = "src",
    doctest = true,
    warnonly = true,
    sitename = "Controlz.jl",
    clean = true,
    pages = ["Controlz" => "index.md",
             "Transfer Functions" => "tfs.md",
             "Simulation" => "sim.md",
             "Visualization" => "viz.md",
             "Control systems" => "controls.md",
             "FAQ" => "faq.md"
             ],
    format = Documenter.HTML(assets = ["assets/flux.css"])
)

deploydocs(
    repo = "github.com/SimonEnsemble/Controlz.jl.git",
    devbranch = "main",
    push_preview=true
)
