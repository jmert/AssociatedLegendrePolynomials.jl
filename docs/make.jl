using Documenter, Legendre

doctest = "--fix"  in ARGS ? :fix :
          "--test" in ARGS ? true : false

DocMeta.setdocmeta!(Legendre, :DocTestSetup, :(using Legendre); recursive=true)

makedocs(
    format = Documenter.HTML(mathengine=MathJax()),
    sitename = "Legendre Polynomials",
    authors = "Justin Willmert",
    modules = [Legendre],
    doctest = doctest,
    doctestfilters = Regex[
        r"Ptr{0x[0-9a-f]+}",
        r"[0-9\.]+ seconds \(.*\)"
    ],
    pages = [
        "Legendre.jl Documentation" => "index.md",
        "Associated Legendre Functions" => [
            "Introduction" => "man/intro.md",
            "Usage" => "man/usage.md",
            "Developer Documentation" => "man/devdocs.md",
            "Literature/References" => "man/references.md"
        ],
        "API Reference" => [
            "Public" => "lib/public.md",
            "Private" => "lib/private.md"
        ]
    ],
    repo = "https://github.com/jmert/Legendre.jl/blob/{commit}{path}#L{line}",
)

deploydocs(
    repo = "github.com/jmert/Legendre.jl.git",
    target = "build",
)
