using BenchmarkTools, Legendre, Markdown

dotune  = "--tune" in ARGS
saveres = "--save" in ARGS

const PARAMS_PATH = joinpath(dirname(@__FILE__), "params.json")
const SUITE   = BenchmarkGroup()
const ASSERTS = BenchmarkGroup()

const LMAX = 700
const MMAX = 350
const NORMS = (LegendreUnitNorm(),
               LegendreSphereNorm(),
               LegendreUnitCoeff{Float64}(LMAX),
               LegendreSphereCoeff{Float64}(LMAX))

const xs = 0.5
const x0 = fill(0.5)
const x1 = collect(range(-1, 1, length=100))
const x2 = collect(reshape(x1, 10, 10))

##############################################################################
# Benchmarks which should be compared revision-to-revision
##############################################################################

@info "Loading benchmarks..."
# Run suite of benchmarks on each normalization style
for norm in NORMS
    nname = string(typeof(norm))
    gnorm = addgroup!(SUITE, nname)

    # Organize by output shape --- single value at (lmax,mmax), vector of values
    # ({(l,mmax) | l <= lmax}), and matrix of values ({(l,m) | l <= lmax and m <= l}).
    for (outdim, osz) in zip(("dim0", "dim1", "dim2"), ((), (LMAX,), (LMAX,MMAX)))
        gname = "out" * outdim
        # Then use true scalar and 0 through 2-dimensional inputs
        for (indim, x) in zip(("scalar", "dim0", "dim1", "dim2"), (xs, x0, x1, x2))
            bname = "in" * indim
            Λ = zeros(size(x)..., (osz .+ 1)...)
            gnorm[(gname, bname)] = @benchmarkable legendre!($norm, $Λ, $LMAX, $MMAX, $x)
        end
    end
end

##############################################################################
# (Optionally tune and) Run the benchmarks
##############################################################################

if dotune
    @info "Tuning benchmarks..."
    tune!(SUITE, verbose=true)
    BenchmarkTools.save(PARAMS_PATH, params(SUITE))
end
if isfile(PARAMS_PATH)
    @info "Loaded benchmarking parameters." PARAMS_PATH
    loadparams!(SUITE, BenchmarkTools.load(PARAMS_PATH)[1], :evals)
end

@info "Running benchmarks..."
res = run(SUITE, verbose=true)
# Truncate down to the minimums immediately. Make the saved file smaller.
res = minimum(res)

##############################################################################
# Assertions on performance characteristics which sould be true for any
# given revision. Uses results of the above benchmarks without comparison to
# another revision.
##############################################################################

@info "Making intra-benchmark comparisons..."
for norm in NORMS
    nname = string(typeof(norm))
    gn = addgroup!(ASSERTS, nname)
    resn = res[nname]

    # Make sure scalar arguments perform better than an equivalent 0-dimensional input
    for outdim in ("dim0", "dim1", "dim2")
        gname = "out" * outdim
        gn[(gname, "inscalar vs indim0")] =
            judge(resn[(gname, "inscalar")], resn[(gname, "indim0")])
    end
end

if saveres
    branch = strip(read(`git symbolic-ref --short HEAD`, String))
    if branch == "master"
        branch = strip(read(`git describe --tags`, String))
    end
    # TODO: BenchmarkTools.load() doesn't work with the judged results, so save the
    #       benchmarks and assertions separately.
    benchfile   = "bench-" * branch * ".json"
    assertsfile = "bench-" * branch * "-asserts.json"
    @info "Saving benchmarks to `bench-$branch.json`" pwd=pwd()
    BenchmarkTools.save(benchfile, res)
    BenchmarkTools.save(assertsfile, ASSERTS)
end

##############################################################################
# Summary printing
##############################################################################

function printtable(group)
    for k in sort!(collect(keys(group)))
        g = group[k]
        md = Markdown.MD()
        push!(md.content, Markdown.parse("## " * k))
        tab = Markdown.Table(Any[], Any[:l, :l])
        push!(tab.rows, Any["Name", "Trial"])
        for l in sort!(leaves(g))
            n, t = l
            push!(tab.rows, Any[n[1], repr(t)])
        end
        push!(md.content, tab)
        show(md)
        print("\n\n")
    end
    print("\n")
end

println("---")
println("# BENCHMARKS")
printtable(res)
println()

println("---")
println("# ASSERTIONS")
printtable(ASSERTS)
println()
