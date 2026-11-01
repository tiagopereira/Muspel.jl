


@testset "zeeman.jl" begin
    @testset "g factors" begin
        @test g_LS(0, 0, 0) == 1
        @test g_LS(1, 1, 0) == 1
        @test g_LS(1, 0, 1) == 2
        @test g_LS(1, 1, 1) == 1.5
        @test g_jj(0, 1, 1, 1, 1) == 1
        @test g_JK(0, 1, 1, 1, 1, 1) == 1
        @test typeof(g_LS(0, 0, 1)) == Float32
        @test typeof(g_LS(0., 0., 1.)) == Float64
        @test typeof(g_jj(0, 0, 1, 1, 1)) == Float32
        @test typeof(g_jj(0., 0., 1., 1., 1.)) == Float64
        @test typeof(g_JK(0, 0, 1, 1, 1, 0)) == Float32
        @test typeof(g_JK(0., 0., 1., 1., 1., 0.)) == Float64
        @test g_eff(0, 0, 1, 1) == 0
        @test g_eff(1, 1, 1, 1) == 1
        @test g_eff(2.5, 1, 2, 1) == g_eff(1, 2.5, 1, 2)
    end
    @testset "zeeman_strength" begin
        @test zeeman_strength(1, 1, -1, 0) ==  zeeman_strength(1, 1, 0, 1)
        @test zeeman_strength(1, 1, 1, 1) == 0.5
        @test zeeman_strength(2, 2, 1, 1) == 0.1
        @test zeeman_strength(1, 0, -1, 0) == zeeman_strength(0, 1, 0, -1)
        # normal zeeman: only one component of each, strength is one
        @test zeeman_strength(1, 0, 1, 0) == 1
        @test zeeman_strength(1, 0, 0, 0) == 1
        @test zeeman_strength(1, 0, -1, 0) == 1
        # two components: sum must add to one
        @test zeeman_strength(1, 1, 1, 0) + zeeman_strength(1, 1, 0, -1) == 1
        @test zeeman_strength(1, 1, 0, 1) + zeeman_strength(1, 1, -1, 0) == 1
        @test zeeman_strength(1, 1, 1, 1) + zeeman_strength(1, 1, -1, -1) == 1
        # normalisation of many components
        tmp_b = tmp_pi = tmp_r = 0
        J = 10
        for Ml in -J:J,  Mu in -J:J
            ΔM = Mu - Ml
            if ΔM == -1
                tmp_r += zeeman_strength(J, J, Ml, Mu)
            elseif ΔM == 0
                tmp_pi += zeeman_strength(J, J, Ml, Mu)
            elseif ΔM == 1
                tmp_b += zeeman_strength(J, J, Ml, Mu)
            end
        end
        @test tmp_b ≈ tmp_pi ≈ tmp_r ≈ 1
        # Test invalid transitions
        @test_throws DomainError zeeman_strength(0, 0, 1, 0)
        @test_throws DomainError zeeman_strength(3, 0, 1, 0)
        @test_throws DomainError zeeman_strength(1, 1, 2, 0)
        @test_throws DomainError zeeman_strength(1, 3, 0, 1)
        @test_throws DomainError zeeman_strength(2, 0, 0, 0)
    end
    @testset "parse_label" begin
        @test Muspel.parse_label_LS("CA II 3P6 4S 2SE") == (0.5, 0)
        @test Muspel.parse_label_LS("CA II 3P6 3D 2DE 3") == (0.5, 2)
        @test Muspel.parse_label_LS("CA II 3P6 3D 2DE 5") == (0.5, 2)
        @test Muspel.parse_label_LS("CA II 3P6 4P 2PO 1") == (0.5, 1)
        @test Muspel.parse_label_LS("CA II 3P6 4P 2PO 3") == (0.5, 1)
        @test Muspel.parse_label_LS("CA III 3P6 1SE") == (0, 0)
        @test Muspel.parse_label_LS("FeI 3d6 4s 5s e5DE 2") == (2, 2)
        @test Muspel.parse_label_LS("FEI 3D6 4S 4P z7DE 3") == (3, 2)
        @test Muspel.parse_label_LS("O II 2P2(3PE) 3P 4PO") == (1.5, 1)
        @test Muspel.parse_label_LS("HE I 1S 5F 3FO 4") == (1, 3)
        @test Muspel.parse_label_LS("9Z") == (4, 20)
        @test_throws ArgumentError Muspel.parse_label_LS("HE I")
        @test_throws ArgumentError Muspel.parse_label_LS("E asdPFS")
    end
    @testset "field_geometry" begin
        @test field_geometry(0f0, 0f0, 0f0) == (0f0, 1f0, 0f0, 1f0, 0f0)
        B, cosγ, sin²γ, cos2χ, sin2χ = field_geometry(0f0, 0f0, -2f0)
        @test (B, cosγ, sin²γ) == (2f0, -1f0, 0f0)
        B, cosγ, sin²γ, cos2χ, sin2χ = field_geometry(3f0, 0f0, 0f0)
        @test (B, cosγ, sin²γ, cos2χ, sin2χ) == (3f0, 0f0, 1f0, 1f0, 0f0)
        B, cosγ, sin²γ, cos2χ, sin2χ = field_geometry(0f0, 3f0, 0f0)
        @test (cos2χ, sin2χ) == (-1f0, 0f0)
        B, cosγ, sin²γ, cos2χ, sin2χ = field_geometry(1f0, 1f0, 0f0)
        @test (cos2χ, sin2χ) == (0f0, 1f0)
        # Compare with explicit angles
        bx, by, bz = 0.3f0, -0.5f0, 0.7f0
        B, cosγ, sin²γ, cos2χ, sin2χ = field_geometry(bx, by, bz)
        γ = acos(bz / sqrt(bx^2 + by^2 + bz^2))
        χ = atan(by, bx)
        @test B ≈ sqrt(bx^2 + by^2 + bz^2)
        @test all((cosγ, sin²γ, cos2χ, sin2χ) .≈ (cos(γ), sin(γ)^2, cos(2χ), sin(2χ)))
        # Field almost along the line of sight is well behaved in single precision
        _, cosγ, sin²γ, _, _ = field_geometry(1f-4, 0f0, 1f0)
        @test sin²γ ≈ 1f-8 rtol=1e-5
    end
    @testset "stokes_profiles" begin
        voigt_itp = create_voigt_itp(LinRange(1f-4, 1f1, 200), LinRange(0f0, 20f0, 400))
        comps = get_zeeman_components(0, 1, 1, 0, 1, 0)
        line = AtomicLine(
            2, 0f0, 0f0, 1, 1, 0f0, 0f0, 0f0, 500f0, 0f0, 1f0, Float32[500, 501], false, true,
            nothing, nothing, LineBroadening(0f0, Muspel.SVector{0, Float32}(), Muspel.SVector{0, Float32}(),
            Muspel.SVector{0, Float32}(), Muspel.SVector{0, Float32}()), comps...
        )
        args = (line, voigt_itp, 0.1f0, 0.5f0, 0.03f0)
        # No field: unpolarised profile and no magneto-optical effects
        ϕI, ϕQ, ϕU, ϕV, ψQ, ψU, ψV = stokes_profiles(args..., 0f0, 0.6f0, 0.3f0, 0.8f0, 0.6f0)
        @test ϕI ≈ real(voigt_itp(0.1f0, 0.5f0)) / (sqrt(π) * 0.03f0)
        @test all(abs.((ϕQ, ϕU, ϕV, ψQ, ψU, ψV)) .< 1e-6 * ϕI)
        # Longitudinal field: V changes sign with the field direction
        ΔλB = 0.005f0
        p1 = stokes_profiles(args..., ΔλB, 1f0, 0f0, 1f0, 0f0)
        p2 = stokes_profiles(args..., ΔλB, -1f0, 0f0, 1f0, 0f0)
        @test p1[4] ≈ -p2[4]
        @test p1[1] ≈ p2[1]
        @test p1[2] ≈ 0 atol=1e-6
        @test p1[3] ≈ 0 atol=1e-6
        # Transverse field: rotating the azimuth by 45° turns Q into U
        q1 = stokes_profiles(args..., ΔλB, 0f0, 1f0, 1f0, 0f0)
        q2 = stokes_profiles(args..., ΔλB, 0f0, 1f0, 0f0, 1f0)
        @test q1[2] ≈ q2[3] rtol=1e-4
        @test abs(q1[3]) < 1e-6
        @test abs(q1[4]) < 1e-6
    end
    @testset "zeeman_components" begin
        # normal Zeeman
        @test get_zeeman_components(0, 1, 1, 0, 1, 0) == ([1], [-1], [1], [0], [1], [1])
        @test get_zeeman_components(0, 1, 0, 0, 1, 1) == ([1], [-1], [1], [0], [1], [1])
        # anomalous Zeeman
        @test get_zeeman_components(1, 1, 1, 1, 1, 1) == (
            [0.5, 0.5], [-1.5, -1.5], [0.5, 0.5], [0, 0], [0.5, 0.5], [1.5, 1.5]
        )
        # consistency checks
        conf1 = [3.5, 3, 4, 2.5, 4, 4]
        conf2 = [9.5, 10, 16, 10.5, 10, 16]
        for conf in [conf1, conf2]
            σr_S, σr_Δ, π_S, π_Δ, σb_S, σb_Δ = get_zeeman_components(conf...)
            @test length(σr_S) == length(σb_S)
            @test maximum(π_Δ) ≈ -minimum(π_Δ)
            @test sum(σr_S) ≈ 1
            @test sum(π_S) ≈ 1
            @test sum(σb_S) ≈ 1
        end
    end

end
