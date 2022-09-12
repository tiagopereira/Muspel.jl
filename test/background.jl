using AtomicData
using Interpolations
using Muspel
using Test
using Unitful

@testset "background.jl" begin
    abund = get_solar_abundances()
    He = read_atom("test_atoms/He_test.yaml")
    H = read_atom("test_atoms/H_test.yaml")
    H_empty = read_atom("test_atoms/H_test_empty.yaml")
    atoms = Vector{AtomicModel}([He, H, H_empty])
    σ_atom_tables = σ_atoms_bf_tables(atoms)
    λ = collect(LinRange(100, 10050, 3))
    nλ = length(λ)
    λ_index = 1:nλ
    log_ne = 15.0:0.05:20.0
    log_temp = 3.0:0.015:6.0
    table = Tables_σ(λ, log_temp, log_ne, atoms)
    table_no_atoms = Tables_σ(λ, log_temp, log_ne, Vector{AtomicModel}())
    temp = [2050.0, 5770.0, 10100.0]
    n = [1.0e15, 7.1e17, 1.0e20]
    @testset "Tables_σ" begin
        # Wrong argument types:
        @test_throws MethodError Tables_σ(λ, collect(log_temp), log_ne, atoms)
        @test_throws MethodError Tables_σ(λ, log_temp, collect(log_ne), atoms)
        @test_throws MethodError Tables_σ(λ, convert(Float32, log_temp), log_ne, atoms)
        @test_throws MethodError Tables_σ(λ, log_temp, convert(Float32, log_ne), atoms)
        @test_throws MethodError Tables_σ([1f0,2f0,3f0], log_temp, log_ne, atoms)
        # Output table shape:
        @test length(table.table_nh) == nλ
        @test length(table.table_ne) == nλ
        # Output table contents against previous implementation:
        @test all([table.table_nh[i](3.5, 17.5) for i in 1:nλ] .≈ [1.2692617030498545e-30,
                                                                  3.5228922025484505e-30,
                                                                  1.2077261039269983e-29])
        @test all([table.table_ne[i](3.5) for i in 1:nλ] .≈ [2.654640736076545e-50,
                                                             2.1537554455941164e-45,
                                                             1.0736587512762791e-44])
    end
    @testset "α_cont_no_atoms" begin
        # Against previous implementation:
        @test all(
            α_cont_no_atoms.(λ, 6e3, 2.0e18, 1.0e20, 1.0e20, 2.0e18) .≈ [
                2.3161098249909214e-10, 6.921500616410073e-9, 2.8583289122838847e-8
            ]
        )
    end
    @testset "α_cont" begin
        # Against previous implementation:
        @test all(α_cont.(
            Ref(Vector{AtomicModel}()),
            Ref(Vector{Vector{Interpolations.FilledExtrapolation}}()),
            λ, 1e4, 2.0e19, 1.0e20, 1.0e20, 2.0e19
        ) .≈ [1.623575604303239e-9, 2.598702182375118e-7, 1.144415527672781e-6])
        @test all(α_cont.(
            Ref(atoms),
            Ref(σ_atom_tables),
            λ, 6e3, 1.0e15, 1.0e20, 1.0e20, 1.0e15
        ) .≈ [7.195575262883325e-12, 1.3271047477761172e-12, 4.3615492198140114e-12])
    end
    @testset "α_cont_fromtables" begin
        # Against previous implementation:
        @test all(α_cont_fromtables.(
                Ref(table), λ_index, λ, 1e4, 2.0e19, 1.0e20, 1.0e20, 2.0e19
            ) .≈ [7.157428857640934e-9, 2.598696308517904e-7, 1.1444077942980045e-6]
        )
        @test all(α_cont_fromtables.(
                Ref(table), λ_index, λ, 6e3, 2.0e15, 1.0e20, 1.0e20, 2.0e15
            ) .≈ [1.2119830339514767e-11, 2.655655166595975e-12, 8.733655789391957e-12]
        )
        @test all(α_cont_fromtables.(
                Ref(table_no_atoms), λ_index, λ, 6e3, 2.0e15, 1.0e20, 1.0e20, 2.0e15
            ) .≈ [2.315341055886709e-13, 2.655655166595975e-12, 8.733655789391957e-12]
        )
        # Test against α_cont and α_cont_no_atoms
        @test all(
            isapprox.(
                α_cont_fromtables.(
                    Ref(table), λ_index, λ, 6e3, 1.0e15,1.0e20, 1.0e20, 1.0e15),
                α_cont.(
                    Ref(atoms), Ref(σ_atom_tables), λ, 6e3, 1.0e15, 1.0e20, 1.0e20, 1.0e15),
                rtol=1e-4
            )
        )
        @test all(
            isapprox.(
                α_cont_fromtables.(Ref(table_no_atoms), λ_index, λ,
                                        6e3, 2.0e18, 1.0e20, 1.0e20, 2.0e18),
                α_cont_no_atoms.(λ, 6e3, 2.0e18, 1.0e20, 1.0e20, 2.0e18),
                rtol=1e-4
            )
        )
    end
    @testset "σ_atoms_bf_tables" begin
        @test length(σ_atom_tables) == length(atoms)
        @test length.(σ_atom_tables) == length.([atom.continua for atom in atoms])
        # Test against values from atomfile * abundances:
        # Using current abundances:
        @test σ_atom_tables[1][1](50.427) / abund[atoms[1].element] == 6.783e-22
        @test σ_atom_tables[1][2](14.609) / abund[atoms[1].element] == 4.844e-23
        @test σ_atom_tables[2][2](364.70521515693) / abund[atoms[2].element] ≈ 1.379e-21
        # Using previous abundances:
        @test σ_atom_tables[1][1](50.427) == 6.783e-22 * 0.08203515443298176
        @test σ_atom_tables[1][2](14.609) == 4.844e-23 * 0.08203515443298176
        @test σ_atom_tables[2][2](364.70521515693) ≈ 1.379e-21 * 1.0
        # A few edges:
        @test σ_atom_tables[1][1](50.428) == 0
        @test σ_atom_tables[1][1](12.7) == 0
        @test σ_atom_tables[2][1](22.793) == 0
        @test σ_atom_tables[2][2](364.8) == 0
        @test σ_atom_tables[2][2](91.175) == 0
    end
    @testset "σ_atoms_bf" begin
        # Against previous implementation:
        @test σ_atoms_bf(σ_atom_tables, atoms, 251.0, 6000.0, 1.0e19) ≈ 5.2185826e-30
        @test σ_atoms_bf(σ_atom_tables, atoms, 10.0, 6000.0, 1.0e19) == 0
        @test σ_atoms_bf(σ_atom_tables, atoms, 370.0, 6000.0, 1.0e19) == 0
    end
    @testset "α_atoms_bf" begin
        # Against previous implementation:
        @test α_atoms_bf(
            σ_atom_tables, atoms, 251.0, 6000.0, 1.0e19, 1.0e20) ≈ 5.2185826e-10
        @test α_atoms_bf(σ_atom_tables, atoms, 10.0, 6000.0, 1.0e19, 1.0e20) == 0
        @test α_atoms_bf(σ_atom_tables, atoms, 370.0, 6000.0, 1.0e19, 1.0e20) == 0
    end
end
