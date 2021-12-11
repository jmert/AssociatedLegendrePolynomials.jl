using SnoopCompile, AssociatedLegendrePolynomials
inf_timing = @snoopi tmin=0.005 begin
    # Invoking aliased wrappers forces inference of underlying implmentation
    for T in (Float64, Float32), I in (Int,)
        s = T(0.5)
        v = collect(range(T(-0.5), T(0.5), length=5))
        Ptab = LegendreUnitCoeff{T}(I(2))
        Λtab = LegendreSphereCoeff{T}(I(2))
        # Scalar argument, single (ℓ,m)
        Plm(I(2), I(2), s)
        Ptab(I(2), I(2), s)
        λlm(I(2), I(2), s)
        Λtab(I(2), I(2), s)
        # Scalar argument, vector ℓ, fixed m
        Plm(I(0):I(2), I(2), s)
        Ptab(I(0):I(2), I(2), s)
        λlm(I(0):I(2), I(2), s)
        Λtab(I(0):I(2), I(2), s)
        # Scalar argument, matrix over all ℓ, m
        Plm(I(0):I(2), I(0):I(2), s)
        Ptab(I(0):I(2), I(0):I(2), s)
        λlm(I(0):I(2), I(0):I(2), s)
        Λtab(I(0):I(2), I(0):I(2), s)

        # Broadcasted operations:
        #  Scalar argument, vector output
        Plm.(I(0):I(2), I(2), s)
        Ptab.(I(0):I(2), I(2), s)
        λlm.(I(0):I(2), I(2), s)
        Λtab.(I(0):I(2), I(2), s)
        #  Scalar argument, matrix output
        Plm.(I(0):I(2), I(0):I(2), s)
        Ptab.(I(0):I(2), I(0):I(2), s)
        λlm.(I(0):I(2), I(0):I(2), s)
        Λtab.(I(0):I(2), I(0):I(2), s)
        #  Vector argument, vector output
        Plm.(I(0):I(2), I(2), v)
        Ptab.(I(0):I(2), I(2), v)
        λlm.(I(0):I(2), I(2), v)
        Λtab.(I(0):I(2), I(2), v)
        #  Vector argument, matrix output
        Plm.(I(0):I(2), I(0):I(2), v)
        Ptab.(I(0):I(2), I(0):I(2), v)
        λlm.(I(0):I(2), I(0):I(2), v)
        Ptab.(I(0):I(2), I(0):I(2), v)

        # Force precompilation of displaying functions
        display(Ptab)
        display(Λtab)
    end
end

pc = SnoopCompile.parcel(inf_timing)
pc = filter!(p -> p.first === :AssociatedLegendrePolynomials, pc)
sort!(pc[:AssociatedLegendrePolynomials])

outfile = joinpath(tempdir(), "precompile")
println("Saving precompile scripts to $outfile")
SnoopCompile.write(outfile, pc)
