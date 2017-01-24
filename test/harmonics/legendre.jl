module Legendre
    using CMB.Harmonics
    using Base.Test

    # In general, all analytically-defined answers are computed using
    # extended-precision ("big") floats and ints. This provides *some* testing
    # for numerical accuracy, but explicit tests of numerical accuracy should
    # probably also be defined.


    #######################
    # LEGENDRE POLYNOMIALS
    #######################

    doc"""
    $P_0$ is a constant. Verify the output is invariant.
    """
    function Pl_00_invariance()
        pl = Harmonics.Pl(0, 2*rand()-1);
        @test pl[1] == 1.0
    end

    ################################
    # ASSOCIATED LEGENDRE FUNCTIONS
    ################################

    doc"""
    $P_0^0$ is constant. Verify the output is invariant for several inputs.
    """
    function Plm_00_invariance()
        pp = Harmonics.PlmPlan(Float64, 1);
        plm = Harmonics.Plm(pp, 2*rand()-1);
        @test plm[1,1] == 1.0
    end

    doc"""
    Match $P_ℓ^{m=0}$ terms for $\ell=1...9$ where the analytical expressions
    have been taken from a table of equations, assuming $x= \cos(\pi/4)$.
    """
    function Plm_l0_analytic()
        pp = Harmonics.PlmPlan(Float64,10)
        plm = Harmonics.Plm(pp, 0.5*sqrt(2));
        @test plm[2,1]  ≈  sqrt(2)/2
        @test plm[3,1]  ≈  1/4
        @test plm[4,1]  ≈ -sqrt(2)/8
        @test plm[5,1]  ≈ -13/32
        @test plm[6,1]  ≈ -17sqrt(2)/64
        @test plm[7,1]  ≈ -19/128
        @test plm[8,1]  ≈  23sqrt(2)/256
        @test plm[9,1]  ≈  611/2048
        @test plm[10,1]  ≈  827sqrt(2)/4096
    end

    doc"""
    The associated Legendre functions $P_ℓ^m$ should equate to the Legendre
    polynomials $P_ℓ$ for $m=0$, so check this is true.
    """
    function Plm_l0_equals_Pl_l()
        const LMAX = 10;
        pl  = Harmonics.Pl(LMAX, cosd(-57.5));
        pp  = Harmonics.PlmPlan(LMAX);
        plm = Harmonics.Plm(pp, cosd(-57.5));
        @test all(pl .≈ plm[:,1])
    end

    ##############################################################
    # SPHERICAL HARMONIC NORMALIZED ASSOCIATED LEGENDRE FUNCTIONS
    ##############################################################

    doc"""
    The initialization condition

    $P_ℓ^ℓ(x) = (-1)^ℓ (2ℓ-1)!! (1-x^2)^(ℓ/2)$

    can be used with extended-precision computations for special values
    of $ℓ$ and $m$. To account for the spherical harmonic normalization,
    the formula can be simplified to

    $\bar{P}_ℓ^ℓ(x) = (-1)^ℓ \frac{1}{\sqrt{2π}} \sqrt{\frac{(2ℓ+1)!!}{(2ℓ)!!} (1-x^2)^(ℓ/2)$

    We'll compute the coefficient separately since it's "easy", and then
    we can use convenient angles that also expand to an exact analytic
    value for $(1-x^2)^(ℓ/2)$.

    Assuming $x = \cos(θ)$, the $(1-x^2)^(ℓ/2)$ term simplifies to an
    expression of $\sin(θ)$:

    $(1-x^2)^(ℓ/2) = (\sin^2 θ)^(ℓ/2) = (\sin θ)^ℓ$

    For $θ = π/4$:

    * even $ℓ$ => $\frac{1}{2^(ℓ/2)}$
    * odd $ℓ$ => $\frac{1}{2^(k/2) \sqrt{2}}$ where $k=\lfloor ℓ/2\rfloor$.
    """
    function Plm_ll_spotchecks()

        function SphericalPll_coeff(T, ℓ)
            factorials = mapreduce(x->sqrt((2x+1)/(2x)), *, 1.0, 1:ℓ);
            sign = (ℓ % 2 == 0) ? 1 : -1;
            return sign * sqrt(1/2π) * factorials;
        end

        pp = Harmonics.PlmSphericalPlan(Float64, 99);
        plm = Harmonics.Plm(pp, cosd(45.0));

        @test plm[3,3]   ≈ SphericalPll_coeff(Float64,2)  / (2^1)
        @test plm[4,4]   ≈ SphericalPll_coeff(Float64,3)  / (2^1 * sqrt(2))
        @test plm[21,21] ≈ SphericalPll_coeff(Float64,20) / (2^10)
        @test plm[22,22] ≈ SphericalPll_coeff(Float64,21) / (2^10 * sqrt(2))
        @test plm[99,99] ≈ SphericalPll_coeff(Float64,98) / (2^49)
        @test plm[100,100] ≈ SphericalPll_coeff(Float64,99) / (2^49 * sqrt(2))
    end

    doc"""
    The normal $P_ℓ^m$ should be equal to the spherical harmonic normalized
    $\bar{P}_ℓ^m$ if we manually normalize them.
    """
    function NormPlm_equals_SphericalPlm()
        const LMAX = 5;
        pp_norm = Harmonics.PlmPlan(LMAX);
        pp_sphr = Harmonics.PlmSphericalPlan(LMAX);
        plm_norm = Harmonics.Plm(pp_norm, cosd(-57.5));
        plm_sphr = Harmonics.Plm(pp_sphr, cosd(-57.5));
        for ℓ=1:LMAX
            for m=0:(ℓ-1)
                factorials = mapreduce(x->1/sqrt(x), *, 1.0, (ℓ-m+1):(ℓ+m));
                norm = sqrt((2ℓ+1)/(2π)) * factorials;
                @test norm*plm_norm[ℓ+1,m+1] ≈ plm_sphr[ℓ+1,m+1]
            end
        end
    end

    function runtests()
        # LEGENDRE POLYNOMIALS
        Pl_00_invariance()
        # ASSOCIATED LEGENDGRE FUNCTIONS
        Plm_00_invariance()
        Plm_l0_equals_Pl_l()
        # SPHERICAL HARMONIC NORMALIZED ASSOCIATED LEGENDRE FUNCTIONS
        Plm_l0_analytic()
        Plm_ll_spotchecks()
        NormPlm_equals_SphericalPlm()
    end
end

