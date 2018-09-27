using Documenter, IGA

makedocs(
    sitename = "IGA.jl",
    format = :html,
    pages = ["Home" => "index.md",
             "API Documentation" => "api.md"]
)
