module ScalarTest

using Test
using Legendre: Scalar

@test Scalar{Float64}() isa Scalar{Float64}
@test Scalar(1.0) isa Scalar{Float64}
@test Scalar{Float64}(1) isa Scalar{Float64}
@test isassigned(Scalar{Float64}())
@test !isassigned(Scalar{BigFloat}())

Z = CartesianIndex()
I = CartesianIndices(())
s = Scalar{Float64}()
s[] = 1
@test s[] == 1
s[Z] = 2
@test s[Z] == 2
s[I] = fill(3)
@test s[I] == fill(3)
@test size(s) == ()
@test ndims(s) == 0
@test length(s) == 1
@test similar(s) isa Scalar{Float64}
@test similar(s, size(s)) isa Scalar{Float64}

@test [x for x in s] == fill(s[])

s .= fill(4)
@test s[] == 4
@test all(s .+ 1 .== 5)
s[Z] .= fill(5)
@test s[] == 5
@test all(s[Z] .+ 1 .== 6)
s[I] .= fill(6)
@test s[I] == fill(6)
@test all(s[I] .+ 1 .== 7)

end # module ScalarTest
