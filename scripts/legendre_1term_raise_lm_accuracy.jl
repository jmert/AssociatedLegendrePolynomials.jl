"""
    legendre_1term_raise_lm_accuracy.jl

Simple tests to explore the numerical accuracy of the 1-term recurrence relation which
boosts (m,m) -> (m+1,m+1) in the associated Legendre polynomials.

## Tests

- `recur0`: Baseline version which is a direct translation of the math.
- `recur1`: Uses FMA instructions to compute (1-x^2). Eliminating intermediate rounding
  is expected to increase the accuracy of the entire recursion.
- `recur2`: Same as `recur1` but also iterates upward in degree/order in pairs. Stepping
  from even to odd ``(m,m) -> (m+1,m+1)`` multiplies by ``y = \\sqrt{1-x^2}`` and the
  coefficient ``μ_m``. The second half then steps from even to even
  ``(m,m) -> (m+2,m+2)`` by multiplying by ``y^2 = 1-x^2`` (without taking the square
  root) and the pair of coefficients ``μ_m`` and ``μ_{m+1}``.
- `recur3`: Same as `recur2` but replaces the FMA calculation of ``1-x^2`` with the
  product ``(1-x)(1+x)``.

## Conclusions
- Without FMA, the product ``(1-x)(1+x)`` is better for ``|x| ≈ 1`` than ``1 - x^2``.
- FMA is the best option for computing ``y = \\sqrt(1-x^2)`` and should be preferred.
- The pair-wise iteration improves numerical accuracy as ``m \\to \\finty``.
"""

using PyPlot, PyCall

function recur0(x, nstep)
    y = sqrt(one(x) - x^2)
    vals = zeros(typeof(x), 2nstep)
    z = inv(sqrt(4oftype(x, π)))
    @inbounds for ii = 1:2nstep
        μ = sqrt(one(x) + one(x)/2(ii+1))
        z = -μ * y * z
        vals[ii] = z
    end
    return vals
end

function recur1(x, nstep)
    y = sqrt(-fma(x, x, -one(x)))
    vals = zeros(typeof(x), 2nstep)
    z = inv(sqrt(4oftype(x, π)))
    @inbounds for ii = 1:2nstep
        μ = sqrt(one(x) + one(x)/2(ii+1))
        z = -μ * y * z
        vals[ii] = z
    end
    return vals
end

function recur2(x, nstep)
    y² = -fma(x, x, -one(x))
    y¹ = sqrt(y²)
    vals = zeros(typeof(x), 2nstep)
    z = inv(sqrt(4oftype(x, π)))
    @inbounds for ii = 1:2:2nstep
        μ₁ = sqrt(one(x) + one(x)/2(ii+1))
        μ₂ = sqrt(one(x) + one(x)/2(ii+2))
        vals[ii] = -μ₁ * y¹ * z
        z = μ₁ * μ₂ * y² * z
        vals[ii+1] = z
    end
    return vals
end

function recur3(x, nstep)
    y² = (one(x) - x) * (one(x) + x)
    y¹ = sqrt(y²)
    vals = zeros(typeof(x), 2nstep)
    z = inv(sqrt(4oftype(x, π)))
    @inbounds for ii = 1:2:2nstep
        μ₁ = sqrt(one(x) + one(x)/2(ii+1))
        μ₂ = sqrt(one(x) + one(x)/2(ii+2))
        vals[ii] = -μ₁ * y¹ * z
        z = μ₁ * μ₂ * y² * z
        vals[ii+1] = z
    end
    return vals
end

# Half of the lmax/mmax to iterate to.
nstep = 512
# Range of arguments, taken more densely from x ≈ 1 where catastrophic cancellation is
# expected without FMA.
x₀ = 1 .- exp2.(range(-16, -eps(), length=2560))

# Reference values calculated with extended precision.
ref = Float64.(reduce(hcat, setprecision(512) do
    recur1.(big.(x₀), nstep)
end))

v0 = reduce(hcat, recur0.(x₀, nstep))
v1 = reduce(hcat, recur1.(x₀, nstep))
v2 = reduce(hcat, recur2.(x₀, nstep))
v3 = reduce(hcat, recur3.(x₀, nstep))

colors = pyimport("matplotlib.colors")
fig, axs = subplots(2, 2, num = "legendre_1term_raise_lm_accuracy", clear = true)
function plot_compare(ax, v; title=nothing, kws...)
    lnorm = colors.LogNorm(vmin = 1, vmax=1e3)
    ulps = abs.(v .- ref) ./ eps.(ref)
    img = ax.pcolormesh(1 .- x₀, 1:2nstep, ulps,
                    norm=lnorm;
                    kws...)
    ax.set_xscale("log")
    ax.set_xlabel("argument \$(1-x)\$")
    ax.set_ylabel("degree/order")
    cbar = ax.figure.colorbar(img, ax=ax)
    cbar.set_label("ulps");
    title !== nothing && ax.set_title(title)
    ax.text(0.05, 0.9, "ulp max = $(round(maximum(ulps), digits=2))", transform=ax.transAxes)
    nothing
end
plot_compare(axs[1,1], v0, title="v0 — Baseline")
plot_compare(axs[1,2], v1, title="v1 — FMA \$(1 - x^2)\$")
plot_compare(axs[2,1], v2, title="v2 — FMA + pairwise iteration")
plot_compare(axs[2,2], v3, title="v3 — \$(1-x)(1+x)\$ + pairwise iteration")

