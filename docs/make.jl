using Documenter, AssociatedLegendrePolynomials

doctest = "--fix"  in ARGS ? :fix :
          "--test" in ARGS ? true : false

DocMeta.setdocmeta!(AssociatedLegendrePolynomials, :DocTestSetup, :(using AssociatedLegendrePolynomials); recursive=true)

makedocs(
    format = Documenter.HTML(mathengine = Documenter.MathJax3()),
    sitename = "Associated Legendre Polynomials",
    authors = "Justin Willmert",
    modules = [AssociatedLegendrePolynomials],
    doctest = doctest,
    doctestfilters = Regex[
        r"Ptr{0x[0-9a-f]+}",
        r"[0-9\.]+ seconds \(.*\)"
    ],
    pages = [
        "AssociatedLegendrePolynomials.jl Documentation" => "index.md",
        "Associated Legendre Functions" => [
            "Introduction" => "man/intro.md",
            "Usage" => "man/usage.md",
            "Developer Documentation" => "man/devdocs.md",
            "Literature/References" => "man/references.md"
        ],
        "API Reference" => "lib/public.md",
    ],
)

deploydocs(
    repo = "github.com/jmert/AssociatedLegendrePolynomials.jl.git",
    push_preview = true
)
