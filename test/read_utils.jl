using Muspel
using Test
using Unitful
using YAML

@testset "read_utils.jl" begin
    using Transparency: h, c_0
    @testset "read_atom" begin
        # Type
        @test (
            typeof(read_atom("test_atoms/H_test_empty.yaml")) == AtomicModel{4, Float64, Int64}
        )
        @test (typeof(read_atom("test_atoms/H_test.yaml")) == AtomicModel{4, Float64, Int64})
        # Fields in output
        He = read_atom("test_atoms/He_test.yaml")
        H = read_atom("test_atoms/H_test.yaml")
        H_empty = read_atom("test_atoms/H_test_empty.yaml")
        @test He.element == :He
        @test H.element == :H
        @test H.nlevels == 4
        @test He.nlevels == 5
        @test H.nlines == 2
        @test H_empty.nlines == 0
        @test H.ncontinua == 2
        @test H_empty.ncontinua == 0
        @test H.Z == 1
        @test He.Z == 2
        @test H.mass == 1.6738233791328e-27
        @test He.mass == 6.646477321159107e-27
        @test length(H.χ) == H.nlevels
        @test length(H_empty.χ) == H_empty.nlevels
        @test length(He.χ) == He.nlevels
        @test H.χ[1] == 0.0
        @test H.χ[2] ≈ 1.6340148841677004e-18
        @test length(H.g) == H.nlevels
        @test length(H_empty.g) == H_empty.nlevels
        @test H.g[1] == 2
        @test H.g[2] == 8
        @test length(H.stage) == H.nlevels
        @test length(H_empty.stage) == H_empty.nlevels
        @test typeof(H.label) == Vector{String}
        # AtomicLine/AtomicContinuum tested in types.jl
        @test typeof(H.lines) == Vector{AtomicLine}
        @test length(H.lines) == H.nlines
        @test typeof(H_empty.lines) == Vector{AtomicLine}
        @test length(H_empty.lines) == H_empty.nlines
        @test typeof(H.continua) == Vector{AtomicContinuum}
        @test length(H.continua) == H.ncontinua
        @test typeof(H_empty.continua) == Vector{AtomicContinuum}
        @test length(H_empty.continua) == H_empty.ncontinua
    end
    @testset "read_continuum" begin
        tr = ["lev1", "lev1_ion1"]
        χ = [10., 10.5, 11.0] * u"aJ"
        stage = [1, 1, 2]
        level_ids = ["lev1", "lev2", "lev1_ion1"]
        # cross_section
        cs = Dict("value"=>[[1e2, 1e-20],
                            [3e2, 3e-20],
                            [2e2, 2e-20]],
                  "unit"=>["nm", "m^2"],
        )
        cont = Dict("cross_section"=>cs, "transition"=> tr)
        continuum = Muspel.read_continuum(cont, χ, stage, level_ids)
        ## Sorted wavelengths and crossections
        @test all((continuum.λ[2:end] .- continuum.λ[1:end-1]) .> 0)
        for i in 1:3
            idx = findall(continuum.λ .== cs["value"][i][1])
            @test continuum.σ[idx] == [cs["value"][i][2]]
        end
        ## Lengths
        @test continuum.nλ == 3
        @test length(continuum.λ) == 3
        @test length(continuum.σ) == 3
        ## λedge from previous implementation
        @test continuum.λedge ≈ 198.64458571489286
        ## indexes
        @test (continuum.up, continuum.lo) == (3, 1)
        # cross_section_hydrogenic
        cs2 = Dict("nλ"=>3,
                  "λ_min"=>Dict("value"=>1980., "unit"=>"Å"),
                  "σ_peak"=>Dict("value"=>1e-20, "unit"=>"m^2"),
        )
        cont2 = Dict("cross_section_hydrogenic"=>cs2, "transition"=> tr)
        continuum = Muspel.read_continuum(cont2, χ, stage, level_ids)
        ## λmin < λedge
        @test_throws AssertionError Muspel.read_continuum(
                                        cont2, [10.0, 11.0, 12.0] * u"aJ", stage, level_ids)
        ## Lengths
        @test continuum.nλ == 3
        @test length(continuum.λ) == 3
        @test length(continuum.σ) == 3
        ## Values from previous implementation
        @test continuum.λedge ≈ 198.64458571489286
        @test all(continuum.λ .≈ [198.0, 198.32229285744643, 198.64458571489286])
        @test all(continuum.σ .≈ [9.911504012980945e-21, 9.955692900228031e-21, 1.0e-20])
        # Missing keyword
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
        # General
        @test line.χup == 10.5
        @test line.χlo == 10.0
        @test line.gup == 8
        @test line.glo == 2
        @test line.Aul == dline["γ_rad"]["value"]
        @test line.Blu ≈ 1.5622986534095843e-7
        @test line.Bul ≈ 3.9057466335239606e-8
        @test line.λ0 ≈ 397.2891714297857
        @test line.f_value == dline["f_value"]
        @test line.label_up == "b"
        @test line.label_lo == "a"
        @test all(line.γ_vdW_const .== [2.3e-15])
        @test all(line.γ_vdW_exp .== [0.0])
        @test line.γ_quad_stark_const == 0
        # when "data" in keys # Using no 4
        @test line.nλ == 3
        @test line.λ == [100.0, 200.0, 300.0]
        @test line.PRD == false
        @test line.Voigt == false
        # Else when type = RH # Use no 3 and no 2
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
        # Else when type = MULTI # Use no 1
        dline = data["radiative_bound_bound"][1]
        line = Muspel.read_line(
                    dline, χ, g, stage, level_ids, label, mass, FloatT=Float64, IntT=Int64)
        @test line.nλ == 5
        @test all(line.λ .≈ [397.28917142978565, 397.337592128717,
                        397.3860246364557, 397.43463712667443, 403.7840001532265])
        @test line.PRD == false
        @test line.Voigt == true
        # Else unrecognized types
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
        ##
    end
    @testset "calc_λline_MULTI" begin
        ##
    end
    @testset "_assign_unit" begin
        ##
    end
    @testset "_read_transition" begin
        # Check a sample
        dE = (ustrip(c_0 |> u"nm/s") * h / u"s") |> u"aJ"
        χ = [0., 10, 50, (10 + ustrip(dE))]*u"aJ"
        level_ids = ["y", "a", "b", "c"]
        data = Dict("transition" => ["a", "c"])
        w, up, lo =  Muspel._read_transition(data, χ, level_ids)
        @test w ≈ 1.0u"nm"
        @test (up, lo) == (4, 2)
        # up > lo
        data = Dict("transition" => ["c", "a"])
        w, up, lo =  Muspel._read_transition(data, χ, level_ids)
        @test (up, lo) == (4, 2)
    end
    @testset "_read_vdW_single" begin
        ##
    end
    @testset "_read_vdW" begin
        ##
    end
    @testset "_read_quadratic_stark" begin
        ##
    end
end
