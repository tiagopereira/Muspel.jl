using Muspel
using Test
using Transparency: const_quadratic_stark, calc_Aji
using Unitful
using YAML
import PhysicalConstants.CODATA2018: h, c_0

@testset "read_utils.jl" begin
    @testset "read_atom" begin
        # Single/Double precision
        @test (typeof(read_atom("test_atoms/H_test.yaml")) == AtomicModel{4, Float64,
                                                                          Int64})
        @test (typeof(read_atom("test_atoms/H_test.yaml", FloatT=Float32, IntT=Int32)
                     ) == AtomicModel{4, Float32, Int32}
        )
        # Simple fields
        H = read_atom("test_atoms/H_test.yaml")
        H_yaml = YAML.load_file("test_atoms/H_test.yaml")
        H_empty = read_atom("test_atoms/H_test_empty.yaml")
        @test H.element == :H
        @test H.nlevels == 4
        @test H.nlines == 2
        @test H_empty.nlines == H_empty.ncontinua == 0
        @test H.ncontinua == 2
        @test H.Z == 1
        @test H.mass == 1.6738233791328e-27
        @test H_empty.mass == 1.674e-27
        @test H.χ[1] == 0.0
        @test all(H.χ[2:end] ≈ [1.6340148841677004e-18, 1.9366102411805722e-18,
                                2.1786865188450865e-18])
        @test H_empty.χ == H.χ
        @test H.g == [2, 8, 18, 1]
        @test H_empty.g == H.g
        @test H.stage == [1, 1, 1, 2]
        @test H_empty.stage == H.stage
        @test H.label == ["H I 1S 2SE", "H I 2P 2PO", "H I 3D 2DE", "H II"]
        @test H_empty.label == H.label
        # Lines
        # Test fields of H.lines one at a time. (AtomicLine contains string fields)
        @test H_empty.lines == Vector{AtomicLine}()
        @test typeof(H.lines) == Vector{AtomicLine}
        ll = Muspel.read_line.(
                H_yaml["radiative_bound_bound"],
                Ref(H.χ*u"J"),
                Ref(H.g),
                Ref(H.stage),
                Ref(["lev1", "lev2", "lev3", "lev1_ion1"]),
                Ref(H.label),
                H.mass*u"kg"
        )
        @test all([
            all(getfield.(H.lines, f) .== getfield.(ll, f)) for f in fieldnames(AtomicLine)
        ])
        # Continua
        @test H_empty.continua == Vector{AtomicContinuum}()
        @test length(H_empty.continua) == H_empty.ncontinua
        @test all(
            H.continua .== Muspel.read_continuum.(H_yaml["radiative_bound_free"],
                                                  Ref(H.χ*u"J"),
                                                  Ref(H.stage),
                                                  Ref(["lev1", "lev2", "lev3", "lev1_ion1"])
            )
        )
    end
    @testset "read_continuum" begin
        tr = ["lev1", "lev1_ion1"]
        χ = [10., 10.5, 11.0] * u"aJ"
        stage = [1, 1, 2]
        level_ids = ["lev1", "lev2", "lev1_ion1"]
        ### cross_section:
        cs = Dict("value"=>[[1e2, 1e-20],
                            [3e2, 3e-20],
                            [2e2, 2e-20]],
                  "unit"=>["nm", "m^2"],
        )
        cont = Dict("cross_section"=>cs, "transition"=> tr)
        continuum = Muspel.read_continuum(cont, χ, stage, level_ids)
        @test continuum.nλ == length(continuum.λ)
        # Sorted wavelengths and crossections
        @test all((continuum.λ[2:end] .- continuum.λ[1:end-1]) .> 0)
        for i in 1:3
            idx = findall(continuum.λ .== cs["value"][i][1])
            @test continuum.σ[idx] == [cs["value"][i][2]]
        end
        # λedge from previous implementation
        @test continuum.λedge ≈ 198.64458571489286
        # indexes
        @test (continuum.up, continuum.lo) == (3, 1)
        ### cross_section_hydrogenic:
        cs2 = Dict("nλ"=>3,
                  "λ_min"=>Dict("value"=>1980., "unit"=>"Å"),
                  "σ_peak"=>Dict("value"=>1e-20, "unit"=>"m^2"),
        )
        cont2 = Dict("cross_section_hydrogenic"=>cs2, "transition"=> tr)
        continuum = Muspel.read_continuum(cont2, χ, stage, level_ids)
        @test continuum.nλ == 3
        # λmin < λedge
        @test_throws AssertionError Muspel.read_continuum(
                                        cont2, [10.0, 11.0, 12.0] * u"aJ", stage, level_ids)
        # Values from previous implementation
        @test continuum.λedge ≈ 198.64458571489286
        @test all(continuum.λ .≈ [198.0, 198.32229285744643, 198.64458571489286])
        @test all(continuum.σ .≈ [9.911504012980945e-21, 9.955692900228031e-21, 1.0e-20])
        ### Missing keyword:
        @test_throws ErrorException Muspel.read_continuum(
                            Dict("transition"=> ["lev1", "lev1_ion1"]), χ, stage, level_ids)
    end
    @testset "read_line" begin
        tr = ["lev1", "lev2"]
        χ = [10., 10.5, 11.0] * u"aJ"
        stage = [1, 1, 2]
        g = [2, 8, 2]
        level_ids = ["lev1", "lev2", "lev1_ion1"]
        label = ["a", "b", "c"]
        mass = 1e-26u"kg"
        data = YAML.load_file("test_atoms/atom_test.yaml")
        dline = data["radiative_bound_bound"][4]
        line = Muspel.read_line(
                    dline, χ, g, stage, level_ids, label, mass, FloatT=Float64, IntT=Int64)
        # General:
        @test line.χup == 10.5
        @test line.χlo == 10.0
        @test line.gup == 8
        @test line.glo == 2
        @test line.Aul == ustrip(calc_Aji(line.λ0 * u"nm", line.glo / line.gup, line.f_value))
        @test line.Blu ≈ 1.4616286808620544e-8
        @test line.Bul ≈ 3.654071702155136e-9
        @test line.λ0 ≈ 397.2891714297857
        @test line.f_value == dline["f_value"]
        @test line.label_up == "b"
        @test line.label_lo == "a"
        @test line.γ_rad == dline["γ_rad"]["value"]
        @test all(line.γ_hydrogen_const .== [2.3e-15])
        @test all(line.γ_hydrogen_exp .== [0.0])
        @test all(line.γ_electrons_const .== [0.0])
        # when "data" in keys:
        @test line.nλ == 3
        @test line.λ == [100.0, 200.0, 300.0]
        @test line.PRD == false
        @test line.Voigt == false
        # Else when type = RH:
        dline = data["radiative_bound_bound"][3]
        line = Muspel.read_line(
                    dline, χ, g, stage, level_ids, label, mass, FloatT=Float64, IntT=Int64)
        @test line.nλ == 4
        @test all(line.λ .≈ [395.29880147724504, 397.286622522245,
                              397.2917203373264, 399.2795413823264])
        @test line.PRD == false
        @test line.Voigt == false
        dline = data["radiative_bound_bound"][2]
        line = Muspel.read_line(
                    dline, χ, g, stage, level_ids, label, mass, FloatT=Float64, IntT=Int64)
        @test line.nλ == 4
        @test all(line.λ .≈ [397.2915988552346, 397.30589962705903,
                        397.46131378793945, 399.2794199002346])
        @test line.PRD == true
        @test line.Voigt == true
        # Else when type = MULTI:
        dline = data["radiative_bound_bound"][1]
        line = Muspel.read_line(
                    dline, χ, g, stage, level_ids, label, mass, FloatT=Float64, IntT=Int64)
        @test line.nλ == 5
        @test all(line.λ .≈ [397.28917142978565, 397.337592128717,
                        397.3860246364557, 397.43463712667443, 403.7840001532265])
        @test line.PRD == false
        @test line.Voigt == true
        # Else unrecognized types:
        dline_u = dline
        dline_u["type_profile"] = "Something"
        @test_throws ErrorException Muspel.read_line(
                     dline, χ, g, stage, level_ids, label, mass, FloatT=Float64, IntT=Int64)
        dline_u["type_profile"] = "Voigt"
        dline_u["wavelengths"]["type"] = "Something"
        @test_throws ErrorException Muspel.read_line(
                     dline, χ, g, stage, level_ids, label, mass, FloatT=Float64, IntT=Int64)
    end
    @testset "calc_λline_RH" begin
        # Wavelength values are against previous implementation
        λ0 = 500.0u"nm"
        nλ = 4
        qcore = 15.0
        qwing = 600.0
        vξ = 2.5u"km/s"
        # Symm/asymm:
        wav = Muspel.calc_λline_RH(λ0, 5, qcore, qwing, vξ) # Does not include λ0
        @test all(
            ustrip(wav) .≈ [497.4950614115939, 499.99679212558004,
                            500.00320787441996, 502.5049385884061]
        )
        @test all(unit.(wav) .== u"nm")
        wav = Muspel.calc_λline_RH(λ0, nλ, qcore, qwing, vξ; asymm=false)
        @test all(
            ustrip(wav) .≈ [500.0030549856672, 500.02105292375967,
                            500.2166461742894, 502.5047856996533]
        )
        @test all(unit.(wav) .== u"nm")
        # qwing <= 2 * qcore:
        wav = Muspel.calc_λline_RH(λ0, nλ, qcore, 20.0, vξ; asymm=false)
        @test all(
            ustrip(wav) .≈ [500.05559401586635, 500.08339102379955,
                            500.1111880317327, 500.1389850396659]
        )
        # Units:
        wav = Muspel.calc_λline_RH(100.0u"Å", nλ, qcore, qwing, vξ; asymm=false)
        @test all(
            ustrip(wav) .≈ [100.00061099713344, 100.00421058475195,
                            100.04332923485788, 100.50095713993066]
        )
        @test all(unit.(wav) .== u"Å")
        # Single precision:
        @test_throws MethodError Muspel.calc_λline_RH(1f0u"Å", nλ, qcore, qwing, vξ)
        wav = Muspel.calc_λline_RH(1f2u"Å", nλ, convert(Float32, qcore),
                                   convert(Float32, qwing), 2.5f3u"m/s"; asymm=false)
        @test eltype(ustrip(wav)) == Float32
        @test all(
            ustrip(wav) .≈ [100.00061099713344, 100.00421058475195,
                            100.04332923485788, 100.50095713993066]
        )
    end
    @testset "calc_λline_MULTI" begin
        # Wavelength values tested against previous implementation:
        λ0 = 500.0u"nm"
        nλ = 4
        q0 = 3.0
        qmax = 600.0
        vξ = 8.0u"km/s"
        # Symm/asymm, (qmax > q0) and (qmax >= 0) and (q0 >= 0):
        wav = Muspel.calc_λline_MULTI(λ0, nλ, q0, qmax, vξ)
        @test all(
            ustrip(wav) .≈ [492.08474595750886, 499.971569413828,
                            500.02843381973275, 508.17405286233975]
        )
        @test all(unit.(wav) .== u"nm")
        wav = Muspel.calc_λline_MULTI(λ0, 5, q0, qmax, vξ) # Includes λ0
        @test all(
            ustrip(wav) .≈ [492.08474595750886, 499.9523432985373, 500.0,
                            500.0476657878395, 508.17405286233975]
        )
        wav = Muspel.calc_λline_MULTI(λ0, nλ, q0, qmax, vξ; asymm=false)
        @test all(
            ustrip(wav) .≈ [500.0, 500.0284296341159, 500.120260419649, 508.17404850397065]
        )
        # (qmax <= q0) linear spacing:
        wav = Muspel.calc_λline_MULTI(λ0, nλ, q0, 1.0, vξ) # qmax > 0
        @test all(
            ustrip(wav) .≈ [499.98665779223063, 499.99555251829133,
                            500.0044475608306, 500.01334291986547]
        )
        dw = wav[2:end] - wav[1:end-1]
        @test all(isapprox.(dw[1], dw, rtol=1e-3))
        # Single precision:
        @test_throws MethodError Muspel.calc_λline_MULTI(λ0, nλ, 1f0, qmax, vξ)
        wav = Muspel.calc_λline_MULTI(5f3u"Å", nλ, 3.0f1, 6.0f2, 8.0f3u"m/s")
        @test eltype(ustrip(wav)) == Float32
        @test all(
            ustrip(wav .|> u"nm") .≈ [492.08474595750886, 499.971569413828,
                                      500.02843381973275, 508.17405286233975]
        )
    end
    @testset "_assign_unit" begin
        @test Muspel._assign_unit(Dict("value"=>1.0, "unit"=>"km / s")) == 1.0u"km/s"
        @test_throws MethodError Muspel._assign_unit(Dict("value"=>"m", "unit"=>1.0))
    end
    @testset "_read_transition" begin
        dE = (ustrip(c_0 |> u"nm/s") * h / u"s") |> u"aJ"
        χ = [0., 10, 50, (10 + ustrip(dE))]*u"aJ"
        level_ids = ["y", "a", "b", "c"]
        data = Dict("transition" => ["a", "c"])
        w, up, lo =  Muspel._read_transition(data, χ, level_ids)
        @test w ≈ 1.0u"nm"
        @test (up, lo) == (4, 2)
        # Output order, up > lo:
        data = Dict("transition" => ["c", "a"])
        w, up, lo =  Muspel._read_transition(data, χ, level_ids)
        @test (up, lo) == (4, 2)
    end
    mass = 1.67e-27u"kg"
    χup = 1.63u"aJ"
    χlo = 0.0u"aJ"
    χ∞ = 2.18u"aJ"
    Z = 1
    @testset "_read_vdW_single" begin
        # Values tested against previous implementation:
        @test_throws AssertionError Muspel._read_vdW_single(Dict(), mass, χup, χlo, χ∞, Z)
        @test_throws ErrorException Muspel._read_vdW_single(
                                                   Dict("type"=>"s"), mass, χup, χlo, χ∞, Z)
        # type = unsold
        data_1a = Dict("type"=>"Unsold")
        data_1b = Dict("type"=>"unsold", "h_coefficient"=>1.)
        data_1c = Dict("type"=>"unsold", "he_coefficient"=>1.)
        data_1d = Dict("type"=>"unsold", "h_coefficient"=>0.4, "he_coefficient"=>0.5)
        @test Muspel._read_vdW_single(data_1a, mass, χup, χlo, χ∞, Z) == (0.0, 0.3)
        @test all(Muspel._read_vdW_single(
            data_1b, mass, χup, χlo, χ∞, Z) .≈ (1.4630219520027063e-15, 0.3))
        @test all(Muspel._read_vdW_single(
            data_1c, mass, χup, χlo, χ∞, Z) .≈ (1.242488840030318e-16, 0.3))
        @test all(Muspel._read_vdW_single(
            data_1d, mass, χup, χlo, χ∞, Z) .≈ (6.473332228025985e-16, 0.3))
        # type = abo
        data_2 = Dict("type"=>"ABO", "α"=>Dict("value"=>2.0), "σ"=>Dict("value"=>3.0))
        @test all(Muspel._read_vdW_single(
            data_2, mass, χup, χlo, χ∞, Z) .≈ (1.043135275447453e-14, -0.5))
        # type = deridder_rensbergen
        data_3 = Dict(
            "type"=>"deridder_rensbergen",
            "h"=>Dict("α"=>Dict("value"=>2.0, "unit"=>"1e-8*cm^3/s"),"β"=>1.2)
        )
        @test all(Muspel._read_vdW_single(
            data_3, mass, χup, χlo, χ∞, Z) .≈ (4.588496830823747e-14, 1.2))
    end
    @testset "_read_vdW" begin
        d = Dict("type"=>"Unsold")
        @test Muspel._read_vdW(
                Dict("broadening_vanderwaals"=>d), mass, χup, χlo, χ∞, Z) == ([0.0], [0.3])
        @test Muspel._read_vdW(
            Dict("broadening_vanderwaals"=>[d,d]), mass, χup, χlo, χ∞, Z) == ([0.0, 0.0],
                                                                               [0.3, 0.3])
        @test Muspel._read_vdW(Dict("something"=>1), mass, χup, χlo, χ∞, Z) == ([], [])
    end
    @testset "_read_stark_quadratic" begin
        data = Dict()
        @test Muspel._read_stark_quadratic(data, mass, χup, χlo, χ∞, Z) == (Float64[], Float64[])
        data["broadening_stark"] = Dict()
        @test Muspel._read_stark_quadratic(data, mass, χup, χlo, χ∞, Z
            ) == ([ustrip(const_quadratic_stark(mass, χup, χlo, χ∞, Z) |> u"m^3 / s")], [1/6])
        data["broadening_stark"]["C_4"] = Dict("unit"=>"m^3/s", "value"=>2.0)
        @test Muspel._read_stark_quadratic(data, mass, χup, χlo, χ∞, Z) == ([2.0], [0.0])
        data["broadening_stark"]["coefficient"] = 2.0
        @test Muspel._read_stark_quadratic(data, mass, χup, χlo, χ∞, Z) == ([4.0], [0.0])
    end
end
